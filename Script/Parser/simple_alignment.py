import sys
import logging
from Bio import pairwise2

STD_INCLUSION = "CGCCAAAATTGATTCCACTGGAGAAATTGAAAAAAGAAGATCCATTACCTCCACCTGCAAACCGCCCGGGTGGTCGTGCAGCTCTTCCACCACACATTCGTTA" \
                "ANNNNNNNNNNTACGTACTTCTGAGTCCAATTACTCTTCNNNNNNNN"
STD_EXCLUSION = "AGGTGGCACCTGTCGATTCGTCTGCGCCTGTTAAAGAAGATCCATTACCTCCACCTGCAAACCGCCCGGGTGGTCGTGCAGCTCTTCCACCACACATTCGTTA" \
                "ANNNNNNNNNNTACGTACTTCTGAGTCCAATTACTCTTCNNNNNNNN"
STD_UNSPLICED = "GCTTCTGCACTTCTATATGCTCAATCGGTTCTGTCTAAGATCCATTACCTCCACCTGCAAACCGCCCGGGTGGTCGTGCAGCTCTTCCACCACACATTCGTTA" \
                "ANNNNNNNNNNTACGTACTTCTGAGTCCAATTACTCTTCNNNNNNNN"
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
        [barcode, UMI, splicing] = extract_info(sequence_R)
        quality = assess_quality(sequence_F, sequence_R, barcode, UMI, splicing)
        if quality:
            print("{0},{1},{2}".format(barcode, splicing, UMI))
        else:
            logging.debug(
                "%d:\n%s\n%s\nbarcode: %s, UMI: %s, splicing:%s"%(i, sequence_F, sequence_R, barcode, UMI, splicing))
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
    ['GAGAGACATG', 'AGAAAAAA', 0]
    """
    isoform_type = classify_isoform(seq_R)
    alignment = pairwise2.align.globalms(seq_R[-len(STD_TEMPLATE):], STD_TEMPLATE, 2, -1, -2, -.5,penalize_end_gaps=(True, True))
    result = pairwise2.format_alignment(*alignment[0])
    results = result.split("\n")
    bar_index = results[2].find("NNNNNNNNNN")
    barcode = results[0][bar_index: bar_index + 10].strip('-')
    UMI_index = results[2].rfind("NNNNNNNN")
    UMI = results[0][UMI_index:UMI_index + 8].strip('-')
    return [barcode, UMI, isoform_type]

def classify_isoform(seq_R):
    """

    :param seq_R: the reverse read (already reverse complimented)
    :return: int 0 -> exclusion; 1 -> inclusion; 2 -> unspliced; 3 -> unknown
    """
    alignment_include = pairwise2.align.globalms(seq_R, STD_INCLUSION, 2, -1, -2, -.5, score_only=True)
    alignment_exclude = pairwise2.align.globalms(seq_R, STD_EXCLUSION, 2, -1, -2, -.5, score_only=True)
    alignment_unspliced = pairwise2.align.globalms(seq_R, STD_UNSPLICED, 2, -1, -2, -.5, score_only=True)
    lst = [alignment_include, alignment_exclude, alignment_unspliced]
    lst.sort(reverse=True)
    logging.info("%d,%d,%d"%(lst[0], lst[1], lst[2]))
    if lst[0] > 200 and lst[0] - lst[1] >= 50:
        if alignment_include == lst[0]:
            return 1
        elif alignment_exclude == lst[0]:
            return 0
        elif alignment_unspliced == lst[0]:
            return 2
        logging.error("impossible situation: maximum is not in the list")
    elif lst[0] >= 100:
        return 3
    else:
        return -1

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

def main(forward_path, reverse_path):
    logging.basicConfig(level=logging.DEBUG, filename="logging_file_alignment.txt", filemode='w')
    logging.info("execution start")
    extract_seqIO(forward_path, reverse_path)
    logging.info("execution finish")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])