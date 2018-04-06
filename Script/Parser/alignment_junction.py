import sys
import logging
from Bio import pairwise2

STD_INCLUSION = "CGCCAAAATTGATTCCACTGGAGAAATTGAAAAAAGAAGATCCATTACCTCCACCTGCAAACCGCCCGGGTGGTCGTGCAGCTCTTCCACCACACATTCGTTA" \
                "ANNNNNNNNNNTACGTACTTCTGAGTCCAATTACTCTTCNNNNNNNN"
STD_EXCLUSION = "AGGTGGCACCTGTCGATTCGTCTGCGCCTGTTAAAGAAGATCCATTACCTCCACCTGCAAACCGCCCGGGTGGTCGTGCAGCTCTTCCACCACACATTCGTTA" \
                "ANNNNNNNNNNTACGTACTTCTGAGTCCAATTACTCTTCNNNNNNNN"
STD_UNSPLICED = "GCTTCTGCACTTCTATATGCTCAATCGGTTCTGTCTAAGATCCATTACCTCCACCTGCAAACCGCCCGGGTGGTCGTGCAGCTCTTCCACCACACATTCGTTA" \
                "ANNNNNNNNNNTACGTACTTCTGAGTCCAATTACTCTTCNNNNNNNN"
STD_CONTROL_INCLUDSION = 'ATGAGACGCCAAAATTGATTCCACTGGAGAAATTGAAAAAAGAAGATCCATTACCTCCACCTGCAAACCGCCCGGGTGGTCGTGCAGCTCTTCC' \
                         'ACCACACATTCGTTAANNNNTACGTACTTCTGAGTCCAATTACTCTTCNNNNNNNN'
STD_CONTROL_EXCLUDSION = 'TCAGTGAGGTGGCACCTGTCGATTCGTCTGCGCCTGTTAAAGAAGATCCATTACCTCCACCTGCAAACCGCCCGGGTGGTCGTGCAGCTCTTCC' \
                         'ACCACACATTCGTTAANNNNTACGTACTTCTGAGTCCAATTACTCTTCNNNNNNNN'
STD_CONTROL_UNSPLICED = 'CAATCGGTTCTGTCTCTTATAACTCAATTCAAATTTTTTCAGAAGATCCATTACCTCCACCTGCAAACCGCCCGGGTGGTCGTGCAGCTCTTCCA' \
                        'CCACACATTCGTTAANNNNTACGTACTTCTGAGTCCAATTACTCTTCNNNNNNNN'
STD_TEMPLATE = "AAACCGCCCGGGTGGTCGTGCAGCTCTTCCACCACACATTCGTTAANNNNNNNNNNTACGTACTTCTGAGTCCAATTACTCTTCNNNNNNNN"
TEST_SEQ = "GAGGTGGCACCTGTCGATTCGTCTGCGCCTGTTAAAGAAGATCCATTACCTCCACCTGCAAACCGCCCGGGTGGTCGTGCAGCTCTTCCACCACACATTCGTTAAGAG" \
           "AGACATGTACGTACTTCTGAGACCAAATACTCTTCAGAAAAAA"


def extract_seqIO(forward_path, reverse_path):
    '''

    :param forward_path: path of the forward read file
    :param reverse_path: path of the reverse read file
    :return: None
    '''
    forward_file = open(forward_path, 'r')
    reverse_file = open(reverse_path, 'r')
    newline_F = forward_file.readline()
    newline_R = reverse_file.readline()
    i = 1
    while newline_F != '' and newline_R != '':
        sequence_F = newline_F
        sequence_R = reverse_compliment(newline_R)
        try:
            [barcode, UMI, splicing, junction] = extract_info(sequence_R)
            quality = assess_quality(sequence_F, sequence_R, barcode, UMI, splicing)
        except:
            quality = False
            logging.error("unknown error at line {0}".format(i))
        if quality:
            print("{0},{1},{2},{3}".format(barcode, splicing, UMI, junction))
        else:
            logging.debug(
                "%d:\n%s\n%s\nbarcode: %s, UMI: %s, splicing:%s, junction: %s"%(i, sequence_F, sequence_R, barcode, UMI, splicing, junction))
        newline_F = forward_file.readline()
        newline_R = reverse_file.readline()
        i += 1
    if newline_F != '' or newline_R != '':
        logging.error("forward and reverse read file has different line number")
    forward_file.close()
    reverse_file.close()
    logging.info("finish extractSeqIO")

def reverse_compliment(seq):
    """
    Return the reverse complimented seq

    :param seq: the reverse read sequence
    :return: the reverse complimented sequence

    >>> reverse_compliment('GAATTC')
    'GAATTC'
    >>> reverse_compliment('aTGCC')
    'GGCAT'
    >>> reverse_compliment('c')
    'G'
    """
    return_strand = ''
    for nt in seq:
        if nt.upper() == 'A':
            return_strand += 'T'
        if nt.upper() == 'G':
            return_strand += 'C'
        if nt.upper() == 'C':
            return_strand += 'G'
        if nt.upper() == 'T':
            return_strand += 'A'
        if nt.upper() == 'N':
            return_strand == 'N'
        assert 'unexpected nt!!'
    return return_strand[::-1]

def extract_info(seq_R):
    """
    Extract barcode, UMI from the sequence and classify which splicing isoform it is.

    Needs to be Implemented!

    :param seq_R: the reverse read (already reverse complimented)
    :return: [barcode(str), UMI(str), splicing_type(int)]
    >>> seq = 'GAGGTGGCACCTGTCGATTCGTCTGCGCCTGTTAAAGAAGATCCATTACCTCCACCTGCAAACCGCCCGGGTGGTCGTGCAGCTCTTCCACCACACATTCGTTAAGAGAGACATGTACGTACTTCTGAGACCAAATACTCTTCAGAAAAAA'
    >>> extract_info(seq)
    ['GAGAGACATG', 'AGAAAAAA', 0, 1]
    >>> seq2 = 'CAATCGGTTCTGTCTCTTATAACTCAATTCAAATTTTTTCAGAAGATCCATTACCTCCACCTGCAAACCGCCCGGGTGGTCGTGCAGCTCTTCCACCACACATTCGTTAAACTGTACGTACTTCTGAGTCCAATTACTCTTCGGGGGGGG'
    >>> extract_info(seq2)
    ['ACTG', 'GGGGGGGG', 2, 1]
    """
    [isoform_type, junction_present] = classify_isoform(seq_R, STD_INCLUSION, STD_EXCLUSION, STD_UNSPLICED)
    alignment = pairwise2.align.globalms(seq_R[-len(STD_TEMPLATE):], STD_TEMPLATE, 2, -1, -2, -.5,penalize_end_gaps=(True, True))
    result = pairwise2.format_alignment(*alignment[0])
    results = result.split("\n")
    bar_index = results[2].find("NNNNNNNNNN")
    barcode = results[0][bar_index: bar_index + 10].strip('-')
    UMI_index = results[2].rfind("NNNNNNNN")
    UMI = results[0][UMI_index:UMI_index + 8].strip('-')
    if barcode == 'ACTG' or barcode == 'GATC':
        [isoform_type, junction_present] = classify_isoform(seq_R, STD_CONTROL_INCLUDSION, STD_CONTROL_EXCLUDSION, STD_CONTROL_UNSPLICED)
    return [barcode, UMI, isoform_type, junction_present]

def classify_isoform(seq_R, include_template, exclude_template, unspliced_template):
    """

    :param seq_R: the reverse read (already reverse complimented)
    :return: int 0 -> exclusion; 1 -> inclusion; 2 -> unspliced; 3 -> unknown
    """
    alignment_include = pairwise2.align.globalms(seq_R, include_template, 2, -1, -2, -.5, score_only=True)
    alignment_exclude = pairwise2.align.globalms(seq_R, exclude_template, 2, -1, -2, -.5, score_only=True)
    alignment_unspliced = pairwise2.align.globalms(seq_R, unspliced_template, 2, -1, -2, -.5, score_only=True)
    lst = [alignment_include, alignment_exclude, alignment_unspliced]
    lst.sort(reverse=True)
    logging.info("{0},{1},{2}".format(lst[0], lst[1], lst[2]))
    if lst[0] > 200 and lst[0] - lst[1] >= 50:
        if alignment_include == lst[0]:
            isoform = 1
            junction = 1 if 'ATTGAAAAAAGAAGATCCAT' in seq_R else 0
            return [isoform, junction]
        elif alignment_exclude == lst[0]:
            isoform = 0
            junction = 1 if 'GCCTGTTAAAGAAGATCCATT' in seq_R else 0
            return [isoform, junction]
        elif alignment_unspliced == lst[0]:
            isoform = 2
            junction = 1 if 'GGTTCTGTCTCTTATAACTC' in seq_R else 0
            return [isoform, junction]
        logging.error("impossible situation: maximum is not in the list")
    elif lst[0] >= 100:
        return [3,0]
    else:
        return [-1,0]

def assess_quality(forward, reverse, barcode, UMI, splicing):
    """
    Assess the quality of the read

    :param seq: a sequence
    :return: bool

    >>> result = assess_quality("ABCD", "BCDA", "CCCCCCCCCC", "BBBBBBBB", -1)
    >>> result
    False
    >>> result = assess_quality("ABCD", "DCBA", "CCC", "DDDDDDDD", 0)
    >>> result
    False
    """
    result = True
    if splicing < 0:
        result = False
        logging.debug("low compatibility")
    if (len(barcode) != 10 and len(barcode) != 4) or '-' in barcode:
        result = False
        logging.debug("incomplete barcode or wrong barcode length %s"%barcode)
    if len(UMI) != 8 or '-' in UMI:
        result = False
        logging.debug("incomplete UMI or wrong UMI length %s"%UMI)
    if not result:
        alignment = pairwise2.align.globalms(reverse[-len(STD_TEMPLATE):], STD_TEMPLATE, 2, -1, -2, -.5)
        alignment_result = pairwise2.format_alignment(*alignment[0])
        logging.debug("\n%s"%alignment_result)
    return result

def main(forward_path, reverse_path, logging_path):
    logging.basicConfig(level=logging.DEBUG, filename=logging_path, filemode='w')
    logging.info("execution start")
    extract_seqIO(forward_path, reverse_path)
    logging.info("execution finish")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2],sys.argv[3])