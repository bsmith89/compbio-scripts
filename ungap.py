#!/usr/bin/env python3
"""Remove dashes and dots from fasta formatted sequences."""

from Bio.SeqIO import parse, write
import sys
import argparse
from copy import copy
from cli import get_default_parser
import logging


def ungap_recs(records):
    """Ungap sequence in records.

    This generator takes an iterable of Bio.SeqRecord objects and produces
    identical records, but with the sequences ungapped.

    """
    for rec in records:
        rec = copy(rec)
        rec.seq = rec.seq.ungap("-")
        rec.seq = rec.seq.ungap(".")
        yield rec

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                parents=[get_default_parser()])
    p.add_argument('in_handles', nargs='*', type=argparse.FileType('r'),
                   metavar="INFILE",
                   default=[sys.stdin])

    args = p.parse_args()

    logger = logging.getLogger(__name__)
    logging.basicConfig(level=args.log_level)
    logger.debug(args)

    for handle in args.in_handles:
        for rec in ungap_recs(parse(handle, args.fmt_infile)):
            write(trans_rec, sys.stdout, args.fmt_outfile)

if __name__ == '__main__':
    main()
