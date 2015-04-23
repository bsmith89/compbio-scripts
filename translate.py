#!/usr/bin/env python3
"""Translate nucleotide sequences to amino acids."""

from Bio.SeqIO import parse, write
import sys
import argparse
from copy import copy
from cli import get_default_parser
from logging import getLogger


def translate_recs(records):
    """Translate sequence from records.

    This geneator takes a iterable of Bio.SeqRecord objects and the same
    records with sequences translated from nucleotides to amino acids.

    """
    for rec in records:
        # TODO: Make sure that this shallow copy of rec
        # is safe.  Is rec.seq a mutable sequence?
        # No, right?
        rec = copy(rec)
        rec.seq = rec.seq.translate()
        yield rec

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                parents=[get_default_parser()])
    p.add_argument('in_handles', nargs='*', type=argparse.FileType('r'),
                   metavar="INFILE",
                   default=[sys.stdin])

    args = p.parse_args()

    l = getLogger(__name__)
    l.setLevel(args.log_level)
    l.info(args)

    for handle in args.in_handles:
        for trans_rec in translate_recs(parse(handle, args.fmt_infile)):
            write(trans_rec, sys.stdout, args.fmt_outfile)

if __name__ == '__main__':
    main()
