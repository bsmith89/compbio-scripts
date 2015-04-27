#!/usr/bin/env python3
"""Align codons based on a corresponding amino-acid sequence.

Not defensive in the slightest:
    Doesn't check that codons and AAs match.
    Doesn't check that sequence labels match (see "match_order.py").

"""

from Bio.SeqIO import parse, write
from Bio.Seq import Seq
import sys
import argparse
import cli
from copy import copy
import logging


logger = logging.getLogger(__name__)


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

def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__,
                                     parents=[cli.get_base_parser(),
                                              cli.get_seq_out_parser(),
                                              cli.get_align_in_parser(),
                                              cli.get_seq_in_parser(),
                                              ])
    args = parser.parse_args(argv)

    return args

def main():
    args = parse_args(sys.argv[1:])
    logging.basicConfig(level=args.log_level)
    logger.debug(args)

    for rec in backalign_recs(parse(args.in_handle, args.fmt_infile),
                              parse(args.align_handle, args.fmt_align)):
        write(rec, args.out_handle, args.fmt_outfile)

if __name__ == '__main__':
    main()
