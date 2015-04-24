#!/usr/bin/env python3
"""Given a file listing sequence IDs, pull those sequences
from a FASTA file.

"""

from Bio.SeqIO import parse, write
import sys
import argparse
import cli
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
                                parents=[cli.get_default_parser(),
                                         cli.get_infile_parser()])
    p.add_argument('list_handle', type=argparse.FileType('r'),
                   metavar="LISTFILE", help=("list of sequence IDs"))

    args = p.parse_args()

    logger = logging.getLogger(__name__)
    logging.basicConfig(level=args.log_level)
    logger.debug(args)

    for rec in get_recs(parse(args.in_handle, args.fmt_infile),
                        get_list(args.list_handle)):
        write(rec, args.out_handle, args.fmt_outfile)

if __name__ == '__main__':
    main()
