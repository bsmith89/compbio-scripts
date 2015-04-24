#!/usr/bin/env python3
"""Reverse complement nucleotide sequences."""

from Bio.SeqIO import parse, write
import sys
import argparse
from copy import copy
import cli
import logging

def revcompl_recs(records):
    """Reverse complement sequence from records.

    This geneator takes a iterable of Bio.SeqRecord objects and the same
    records with sequences reverse transcribed.

    """
    for rec in records:
        # TODO: Make sure that this shallow copy of rec
        # is safe.  Is rec.seq a mutable sequence?
        # No, right?
        rec = copy(rec)
        rec.seq = rec.seq.reverse_complement()
        yield rec

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                parents=[cli.get_default_parser(),
                                         cli.get_infile_parser()])
    args = p.parse_args()

    logger = logging.getLogger(__name__)
    logging.basicConfig(level=args.log_level)
    logger.debug(args)

    for rc_rec in revcompl_recs(parse(args.in_handle, args.fmt_infile)):
        write(rc_rec, args.out_handle, args.fmt_outfile)

if __name__ == '__main__':
    main()
