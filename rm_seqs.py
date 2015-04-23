#!/usr/bin/env python3
"""Given a list of sequence IDs and a sequence file, output the latter
without any of the former.

"""

from Bio.SeqIO import parse, write
import sys
import argparse
from cli import get_default_parser
from logging import getLogger


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
                                parents=[get_default_parser()])
    p.add_argument('rm_handle', type=argparse.FileType('r'),
                   metavar="LISTFILE")
    p.add_argument('in_handle', nargs='?', type=argparse.FileType('r'),
                   metavar="SEQFILE", default=sys.stdin)

    args = p.parse_args()

    l = getLogger(__name__)
    l.setLevel(args.log_level)
    l.info(args)

    for rec in rm_recs(parse(args.in_handle, args.fmt_infile),
                       get_list(args.rm_handle)):
        write(rec, sys.stdout, args.fmt_outfile)

if __name__ == '__main__':
    main()
