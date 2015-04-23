#!/usr/bin/env python3
"""Given a file listing sequence IDs, pull those IDs from a FASTA file.

"""

from Bio.SeqIO import parse, write
import sys

def main():
    with open(sys.argv[1]) as handle:
        ids = [line.strip() for line in handle]
    for rec in parse(sys.argv[2], 'fasta'):
        if rec.id in ids:
            write(rec, sys.stdout, 'fasta')

if __name__ == '__main__':
    main()
