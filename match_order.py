#!/usr/bin/env python3
"""Given two FASTA files, reorder the second to match IDs with the first.

"""

from Bio.SeqIO import parse, index, write
import sys
import argparse

DEFAULT_INFILE_FMT = 'fasta'
DEFAULT_OUTFILE_FMT = 'fasta'

def order_recs(ordered_recs, recs_map):
    for rec in ordered_recs:
        yield recs_map[rec.id]

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
    p.add_argument('match_path', type=str,
                   metavar="UNORDERED")
    p.add_argument('ord_handle', nargs='?', type=argparse.FileType('r'),
                   metavar="ORDERED", default=sys.stdin)

    args = p.parse_args()

    for rec in order_recs(parse(args.ord_handle, args.fmt_infile),
                          index(args.match_path, args.fmt_infile)):
        write(rec, sys.stdout, args.fmt_outfile)

if __name__ == '__main__':
    main()
