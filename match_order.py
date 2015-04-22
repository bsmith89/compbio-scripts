#!/usr/bin/env python3
"""Given two FASTA files, reorder the second to match IDs with the first.

Outputs to stdout.

"""

from Bio.SeqIO import parse, index, write
import sys

def main():
    ordered_seqs = parse(sys.argv[1], 'fasta')
    unordered_seqs = index(sys.argv[2], 'fasta')

    for seq in ordered_seqs:
        write(unordered_seqs[seq.id], sys.stdout, 'fasta')

if __name__ == '__main__':
    main()
