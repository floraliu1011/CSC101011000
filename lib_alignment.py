from Bio import pairwise2

SEQUENCE = "ACGCGTATGGATTACAAGGATGACGATGACAAGGGGGTACCTGCCCCAAAAAAAAAACGCAAAGTGGAGGAC" \
           "CCAGTACCCGGATCTGAATTCATATATGCTAGCAGTATCAGTGAGGTGGCACCTGTCGATTCGTCTGCGCCT" \
           "GTTAAAGGTAAAAATGTTTGAAAAATTTCTATAGATCAAATATTAAAACCACCCGGACGCGAGATTTTTCCT" \
           "TTTTTTTTGTCCAAAAATCGGTCTCGACACGACAATTTTCGTTATATGCAAACGGATGTGCACCTTTAAAGA" \
           "GTACTGTAAATTAACAAATCGTGAGATAAACTATGAGAAATCGATGAAAATTCCACATCAATGAAACTTTTA" \
           "TATTACAGTACTCTTTAAAGCCACACGCCGAAAAAATTGTAGCGTCGAGACCGGATACCGTATAGACCAAAA" \
           "ATCGCAAAATTTCGCGCCTGGCTAATAAAACACGTTTTTCATTGAGCTGAAATTTGAAAAATTTCGGAAATT" \
           "TATTAGCAATAAATTTTTTAGCTAAAACTAAACTCTCGCGAGTTTCGGAAAATGTTGTGCCGGATGAGACGC" \
           "CAAAATTGATTCCACTGGAGAAATTGAAAAAAGGTGATCCGCAGAGACATGGTTTCAAATTGGNNNNNNNNN" \
           "NTGCTATGTCGTTTTTCAAACCTCTTTACAAACGAAAAATCTCACTATTTGGTTCTAGTCGTATAATTCTTA" \
           "GTTTTATGATTTTCATCCACAAACTTTCAATTTAATTCAAATTTCCATAACTGTCCAATTTTCTTTTAAAAA" \
           "ACTGATTTTAAATGATTCAGTTTCTGCTTCTGCACTTCTATATGCTCAATCGGTTCTGTCTCTTATAACTCA" \
           "ATTCAAATTTTTTCAGAAGATCCATTACCTCCACCTGCAAACCGCCCGGGTGGTCGTGCAGCTCTTCCACCA" \
           "CACATTCGTTAANNNNNNNNNNTACGTA"

def my_alignment(seq, seq2 = SEQUENCE, sliced = True):
    alignment = pairwise2.align.globalms(seq, seq2, 2, -1, -2, -.5)
    result = pairwise2.format_alignment(*alignment[0])
    if sliced:
        results = result.split("\n")
        cis_index = results[3].find("NNNNNNNNNN")
        bar_index = results[3].rfind("NNNNNNNNNN")
        cis = results[0][cis_index:cis_index + 10]
        bar = results[0][bar_index: bar_index + 10]
        return (cis, bar)
    else:
        return result

def alignIO(in_filename, out_filename, sliced=True):
    file = open(in_filename, "r")
    out_file = open(out_filename, "a")
    file.readline()
    seq = file.readline()
    alignment_result = my_alignment(seq, sliced=sliced)
    if type(alignment_result) == tuple:
        index = in_filename[:in_filename.find("-")]
        out_file.write("{0},{1},{2}\n".format(index, alignment_result[0], alignment_result[1]))
        out_file.close()
        file.close()
    else:
        print(alignment_result)

def batch(start, end):
    for i in range(start, end+1):
        in_filename =  str(i) + "-user_added.fasta"
        alignIO(in_filename, "../cis_bar_library.csv")

if __name__ == "__main__":
    import os
    os.chdir("/Users/liuflora/Desktop/ROP299/phase2/library_sanger/")
    batch(11,20)
