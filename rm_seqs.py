#!/usr/bin/env python3
"""Given a list of sequence IDs and a sequence file, output the latter
without any of the former.

File should be FASTA formatted.

"""

from Bio.SeqIO import parse, write
import sys


def main():
    with open(sys.argv[1]) as handle:
        ids = [line.strip() for line in handle]
    for rec in parse(sys.argv[2], 'fasta'):
        if rec.id not in ids:
            write(rec, sys.stdout, 'fasta')

if __name__ == '__main__':
    main()
