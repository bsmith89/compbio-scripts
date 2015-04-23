#!/usr/bin/env python3
"""Given a file listing sequence IDs, pull those sequences
from a FASTA file.

"""

from Bio.SeqIO import parse, write
import sys
import argparse
from cli import get_default_parser
import logging

def get_recs(recs, get_ids):
    for rec in recs:
        if rec.id in get_ids:
            yield rec
        else:
            continue

def get_list(handle):
    out = []
    for line in handle:
        out.append(line.strip())
    return out

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                parents=[get_default_parser()])
    p.add_argument('list_handle', type=argparse.FileType('r'),
                   metavar="LISTFILE")
    p.add_argument('in_handle', nargs='?', type=argparse.FileType('r'),
                   metavar="SEQFILE", default=sys.stdin)

    args = p.parse_args()
    logging.basicConfig(level=args.log_level)
    logging.debug(args)

    for rec in get_recs(parse(args.in_handle, args.fmt_infile),
                       get_list(args.list_handle)):
        write(rec, sys.stdout, args.fmt_outfile)

if __name__ == '__main__':
    main()
