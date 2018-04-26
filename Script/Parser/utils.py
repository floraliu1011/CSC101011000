import logging

def parse(forward_seq, reverse_seq):
    """
    Return the parsed elements if the sequence is qualified, and None if not

    :param forward_seq: the raw sequence in the forward read
    :param reverse_seq: the raw sequence in the reverse read
    :return: a string of the parsed elements
    """
    sequence_F = forward_seq
    sequence_R = reverse_compliment(reverse_seq)
    [barcode, UMI, splicing] = extract_info(sequence_F, sequence_R)
    quality = assess_quality(sequence_F, sequence_R, barcode, UMI, splicing)
    if quality:
        return("{0},{1},{2}".format(barcode, splicing, UMI))
    else:
        logging.debug(
            "%s\n%s\nbarcode: %s, UMI: %s, splicing:%s" % (sequence_F, sequence_R, barcode, UMI, splicing))


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


def extract_info(seq_F, seq_R):
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

    :param forward: the forward read
    :param reverse: the reverse read
    :param barcode: the barcode parsed out
    :param UMI: UMI
    :param splicing: the splicing pattern
    :return: True if the parsed result is qualified
    """
    if (len(barcode) == 10 or len(barcode) == 4) and len(UMI) == 8:
        return True
    return False
