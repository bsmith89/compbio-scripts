#!/usr/bin/env python3
"""Translate nucleotide sequences to amino acids."""

from Bio.SeqIO import parse, write
import sys
import argparse
from copy import copy

DEFAULT_INFILE_FMT = 'fasta'
DEFAULT_OUTFILE_FMT = 'fasta'

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
    p = argparse.ArgumentParser(description=__doc__)
    # Arguments which may be generalizable to many scripts.
    # TODO: Break these out into a module which can be imported.
    p.add_argument("--in-fmt", "-f", dest='fmt_infile', nargs=1, type=str,
                   metavar="FORMAT", default=DEFAULT_INFILE_FMT,
                   help="file format of infile or stdin")
    p.add_argument("--out-fmt", "-t", dest='fmt_outfile', nargs=1, type=str,
                   metavar="FORMAT", default=DEFAULT_OUTFILE_FMT,
                   help="file format of outfile or stdout")
    p.add_argument("--verbose", "-v", action='count',
                   help="increase verbosity")

    # Arguments specific to this script.
    p.add_argument('in_handles', nargs='*', type=argparse.FileType('r'),
                   metavar="INFILE",
                   default=[sys.stdin])

    args = p.parse_args()

    for handle in args.in_handles:
        for trans_rec in translate_recs(parse(handle, args.fmt_infile)):
            write(trans_rec, sys.stdout, args.fmt_outfile)

if __name__ == '__main__':
    main()
