#!/usr/bin/env python3
"""Align codons based on a corresponding amino-acid sequence.

"""

import sys
import argparse
import logging
from warnings import warn

from Bio.SeqIO import parse
from Bio.AlignIO import read as read_aln
from Bio.AlignIO import write
from Bio.Alphabet import generic_protein, generic_nucleotide
from Bio.codonalign import build as build_codon_alignment

import lib.cli as cli

logger = logging.getLogger(__name__)
logging.captureWarnings(True)

def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__,
                                     parents=[cli.get_base_parser(),
                                              cli.get_seq_out_parser(),
                                              cli.get_align_in_parser(),
                                              cli.get_seq_in_parser(),
                                              ])
    args = parser.parse_args(argv[1:])

    return args

def main():
    args = parse_args(sys.argv)
    logging.basicConfig(level=args.log_level)
    logger.debug(args)

    prot_aln = read_aln(args.align_handle, args.fmt_align,
                        alphabet=generic_protein)
    nucl_recs = parse(args.in_handle, args.fmt_infile,
                      alphabet=generic_nucleotide)
    corr_dict = {rec.id: rec.id for rec in prot_aln}

    codon_aln = build_codon_alignment(prot_aln, nucl_recs, corr_dict)
    write(codon_aln, args.out_handle, args.fmt_outfile)

if __name__ == '__main__':
    main()
