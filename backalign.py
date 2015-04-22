#!/usr/bin/env python3
"""Align codons based on a corresponding amino-acid sequence.

Not defensive in the slightest:
    Doesn't check that codons and AAs match
    Doesn't check that sequence labels match (see "match_order.py")

"""

from Bio.SeqIO import parse, write
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys


def codons(sequence):
    """Return groups of 3 items from a sequence.
    
    This is specifically meant to return codons from an in-frame sequence."""
    codon = ""
    for nucl in sequence:
        codon += nucl
        if len(codon) == 3:
            yield codon
            codon = ""

def main():
    for afa_rec, fn_rec in zip(parse(sys.argv[1], 'fasta'),
                               parse(sys.argv[2], 'fasta')):
        afn_str = ""
        codon_iter = codons(fn_rec.seq)
        for amino in afa_rec.seq:
            if amino == ".":
                continue
            elif amino.islower():
                next(codon_iter)
            elif amino == "-":
                afn_str += "---"
            else:
                afn_str += next(codon_iter)
        out_rec = SeqRecord(Seq(afn_str), id=afa_rec.id, description="")
        write(out_rec, sys.stdout, 'fasta')

if __name__ == '__main__':
    main()
