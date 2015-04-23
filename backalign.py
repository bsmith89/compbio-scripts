#!/usr/bin/env python3
"""Align codons based on a corresponding amino-acid sequence.

Not defensive in the slightest:
    Doesn't check that codons and AAs match
    Doesn't check that sequence labels match (see "match_order.py")

"""

from Bio.SeqIO import parse, write
from Bio.Seq import Seq
import sys
import argparse
from cli import get_default_parser
from copy import copy
import logging

def codons(sequence):
    """Return groups of 3 items from a sequence.
    
    This is specifically meant to return codons from an in-frame sequence."""
    codon = ""
    for nucl in sequence:
        codon += nucl
        if len(codon) == 3:
            yield codon
            codon = ""

def backalign(nucl, prot):
    align_nucl = ""
    codon_iter = codons(nucl)
    for amino in prot:
        if amino == ".":
            continue
        elif amino.islower():
            next(codon_iter)
        elif amino == "-":
            align_nucl += "---"
        else:
            align_nucl += next(codon_iter)
    return Seq(align_nucl)

def backalign_recs(nucl_recs, prot_recs):
    for nucl_rec, prot_rec in zip(nucl_recs, prot_recs):
        out_rec = copy(nucl_rec)
        out_rec.seq = backalign(nucl_rec.seq, prot_rec.seq)
        yield out_rec

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                parents=[get_default_parser()])
    p.add_argument('in_prot', type=argparse.FileType('r'),
                   metavar="PROT-FILE")
    p.add_argument('in_nucl', nargs='?', type=argparse.FileType('r'),
                   metavar="NUCL-FILE",
                   default=sys.stdin)

    args = p.parse_args()
    logging.basicConfig(level=args.log_level)
    logging.debug(args)

    for rec in backalign_recs(parse(args.in_nucl, args.fmt_infile),
                              parse(args.in_prot, args.fmt_infile)):
        write(rec, sys.stdout, args.fmt_outfile)

if __name__ == '__main__':
    main()
