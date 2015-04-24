#!/usr/bin/env python3
"""Given two FASTA files, reorder the second to match IDs with the first.

"""

from Bio.SeqIO import parse, index, write
import sys
import argparse
import cli
import logging


def order_recs(ordered_recs, recs_map):
    for rec in ordered_recs:
        yield recs_map[rec.id]

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                parents=[cli.get_default_parser()])
    p.add_argument('match_path', type=str,
                   metavar="UNORDERED", help=("sequences to be put in order"))
    p.add_argument('ord_handle', nargs='?', type=argparse.FileType('r'),
                   metavar="ORDERED", default=sys.stdin,
                   help=("sequences already in order"))

    args = p.parse_args()

    logger = logging.getLogger(__name__)
    logging.basicConfig(level=args.log_level)
    logger.debug(args)

    for rec in order_recs(parse(args.ord_handle, args.fmt_infile),
                          index(args.match_path, args.fmt_infile)):
        write(rec, args.out_handle, args.fmt_outfile)

if __name__ == '__main__':
    main()
