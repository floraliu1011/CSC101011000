import sys
import logging


def extract_seqIO(forward_path, reverse_path):
    '''

    :param forward_path: path of the forward read file
    :param reverse_path: path of the reverse read file
    :return: None
    '''
    logging.info("start extractSeqIO")
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
            logging.debug("%d:\n%s\n%s\nbarcode: %s, UMI: %s, splicing:%s"%(i, sequence_F, sequence_R, barcode, UMI, splicing))
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
    UMI = seq_R[-8:]
    isoform_type = classify_isoform(seq_R)
    barcode = seq_R[seq_R.find('ACATTCGTTAA') + len('ACATTCGTTAA'):seq_R.rfind('TACGTACT')]
    return [barcode, UMI, isoform_type]


def classify_isoform(seq_R):
    """

    :param seq_R: the reverse read (already reverse complimented)
    :return: int 0 -> exclusion; 1 -> inclusion; 2 -> unspliced; 3 -> unknown
    """
    if 'ATTGAAAAAAGAAGATCCAT' in seq_R:
        return 1
    elif 'GCCTGTTAAAGAAGATCCATT' in seq_R:
        return 0
    elif 'GGTTCTGTCTCTTATAACTC' in seq_R:
        return 2
    else:
        return 3


def assess_quality(forward, reverse, barcode, UMI, splicing):
    """
    Assess the quality of the read

    :param seq: a sequence
    :return: bool
    """
    if (len(barcode) == 10 or len(barcode) == 4) and len(UMI) == 8:
        return True
    return False


def main(forward_path, reverse_path):
    logging.basicConfig(level=logging.DEBUG, filename="logging_file.txt", filemode='w')
    logging.info("execution start")
    extract_seqIO(forward_path, reverse_path)
    logging.info("execution finish")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])