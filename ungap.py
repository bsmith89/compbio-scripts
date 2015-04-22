#!/usr/bin/env python3
"""Remove dashes and dots from fasta formatted sequences."""

from Bio.SeqIO import parse, write
import sys

for rec in parse(sys.stdin, 'fasta'):
    rec.seq = rec.seq.ungap("-")
    rec.seq = rec.seq.ungap(".")
    write(rec, sys.stdout, 'fasta')
