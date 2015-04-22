#!/usr/bin/env python3
"""Translate nucleotide sequences to amino-acids.

"""

from Bio.SeqIO import parse, write
import sys

def main():
    for rec in parse(sys.stdin, 'fasta'):
        rec.seq = rec.seq.translate()
        write(rec, sys.stdout, 'fasta')

if __name__ == '__main__':
    main()
