#!/usr/bin/env python3
"""Given a file listing sequence IDs, pull those sequences
from a FASTA file.

"""

from Bio.SeqIO import parse, write
import sys
import argparse
from ..lib import cli
import logging

logger = logging.getLogger(__name__)
logging.captureWarnings(True)


def fetch_iter(rec_iter, id_set):
    for rec in rec_iter:
        if id_set and rec.id in id_set:
            id_set -= {rec.id}
            yield rec

def exclude_iter(rec_iter, id_set):
    return (rec for rec in rec_iter if rec.id not in id_set)

def order_iter(rec_iter, order):
    index = {rec.id: rec for rec in rec_iter}
    return (index[rec_id] for rec_id in order)

def _get_extra_args():
    p = argparse.ArgumentParser(add_help=False)
    g = p.add_argument_group(*cli.FMT_GROUP)
    g.add_argument("-m", "--match-order", action='store_true',
                   help="output in the same order as LIST")
    g.add_argument("-x", "--excluding", action='store_true',
                   help="only output sequences *not* in LIST")
    return p

def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__,
                                     parents=[cli.get_base_parser(),
                                              cli.get_seq_out_parser(),
                                              cli.get_list_in_parser(),
                                              cli.get_seq_in_parser(),
                                              _get_extra_args(),
                                              ])
    args = parser.parse_args(argv[1:])
    return args

def main():
    args = parse_args(sys.argv)
    logging.basicConfig(level=args.log_level)
    logger.debug(args)

    if args.match_order and args.excluding:
        raise ValueError("--match-order and --excluding cannot both be set.")

    to_fetch = [line.strip() for line in args.list_handle]
    fetch_set = set(to_fetch)

    rec_iter = parse(args.in_handle, args.fmt_infile)

    if args.excluding:
        out_iter = exclude_iter(rec_iter, fetch_set)
    elif args.match_order:
        out_iter = order_iter(fetch_iter(rec_iter, fetch_set), to_fetch)
    else:
        out_iter = fetch_iter(rec_iter, fetch_set)

    write(out_iter, sys.stdout, args.fmt_outfile)

if __name__ == '__main__':
    main()
