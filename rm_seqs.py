#!/usr/bin/env python3
"""Given a list of sequence IDs and a sequence file, output the latter
without any of the former.

"""

from Bio.SeqIO import parse, write
import sys
import argparse
import cli
import logging


def rm_recs(recs, rm_ids):
    for rec in recs:
        if rec.id in rm_ids:
            continue
        else:
            yield rec

def get_list(handle):
    out = []
    for line in handle:
        out.append(line.strip())
    return out

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                parents=[cli.get_default_parser(),
                                         cli.get_infile_parser()])
    p.add_argument('rm_handle', type=argparse.FileType('r'),
                   metavar="LISTFILE", help=("list of sequences"))

    args = p.parse_args()

    logger = logging.getLogger(__name__)
    logging.basicConfig(level=args.log_level)
    logger.debug(args)

    for rec in rm_recs(parse(args.in_handle, args.fmt_infile),
                       get_list(args.rm_handle)):
        write(rec, args.out_handle, args.fmt_outfile)

if __name__ == '__main__':
    main()
