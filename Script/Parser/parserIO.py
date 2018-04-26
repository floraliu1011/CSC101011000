import sys
import logging
import utils
import argparse

def parse_seqIO(forward_path, reverse_path, fastq):
    forward_file = open(forward_path, 'r')
    reverse_file = open(reverse_path, 'r')
    if fastq:
        forward_file.readline()
        reverse_file.readline()
    newline_F = forward_file.readline()
    newline_R = reverse_file.readline()
    while newline_F != '' and newline_R != '':
        result = utils.parse(newline_F, newline_R)
        if result:
            print(result)
        if fastq:
            forward_file.readline()
            reverse_file.readline()
            forward_file.readline()
            reverse_file.readline()
            forward_file.readline()
            reverse_file.readline()
        newline_F = forward_file.readline()
        newline_R = reverse_file.readline()
    if newline_F != '' or newline_R != '':
        logging.error("forward and reverse read file has different line number")
    forward_file.close()
    reverse_file.close()
    logging.info("finish extractSeqIO")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A fastq sequence parser.")
    parser.add_argument("-s", action='store_false',
                        help='all lines in the files are DNA sequences')
    parser.add_argument("forward_path", type=str,
                        help="the relative path to the forward read file")
    parser.add_argument("reverse_path", type=str,
                        help="the relative path to the reverse read file")
    parser.add_argument("logging_path", nargs='?', type=str,
                        default='logging_file.txt',
                        help="path to the logging file")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG, filename=args.logging_path, filemode='w')
    logging.info("execution start")
    parse_seqIO(args.forward_path, args.reverse_path, args.s)
    print(args.s)
    logging.info("execution finish")