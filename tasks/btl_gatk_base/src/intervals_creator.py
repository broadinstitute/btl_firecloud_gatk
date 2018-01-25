#!/usr/bin/env python
__author__ = "Amr Abouelleil"

import sys
import argparse
from math import ceil
from Bio import SeqIO

parser = argparse.ArgumentParser(description="")
parser.add_argument('-r', '--reference', action='store', help='Input reference file to create intervals for.')
parser.add_argument('-i', '--interval_size', action='store', help='Maximum size for each interval.', default=1.5)
args = vars(parser.parse_args())


def main():
    interval_bases = args['interval_size'] * 1000000
    for record in SeqIO.parse(args['reference'], 'fasta'):
        count = 0
        if len(record) <= interval_bases:
            print("%s:1-%i" % (record.id, len(record)))
        else:
            chunks = int(ceil(len(record)/interval_bases))
            chunk_start = 1
            chunk_end = interval_bases
            while count < chunks:
                print("%s:%i-%i" % (record.id, chunk_start, chunk_end))
                chunk_start = chunk_end + 1
                chunk_end += interval_bases
                if chunk_end > len(record):
                    chunk_end = len(record)
                count += 1


if __name__ == "__main__":
    sys.exit(main())

