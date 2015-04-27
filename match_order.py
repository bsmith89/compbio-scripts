#!/usr/bin/env python3
"""Given two FASTA files, reorder the second to match IDs with the first.

"""
# TODO: The logic of this script is off.
# Take sequences from STDIN and order them based on a list of sequence IDs
# I could re-write `get_seqs.py` to have an option to respect order and
# then write another script to pull the IDs from a sequence file in order.

from Bio.SeqIO import parse, index, write
import sys
import argparse
import cli
import logging

logger = logging.getLogger(__name__)

def order_recs(ordered_recs, recs_map):
    for rec in ordered_recs:
        yield recs_map[rec.id]

def match_order_parser():
    p = argparse.ArgumentParser(add_help=False)
    g = p.add_argument_group(*cli.POS_GROUP)
    g.add_argument("match_path", type=str, metavar="SEQUENCES",
                   help="sequence file")
    g.add_argument("ord_handle", nargs='?', type=argparse.FileType('r'),
                   metavar="ORDER", default=sys.stdin,
                   help="ordered sequences")
    h = p.add_argument_group(*cli)
    g.add_argument("-f", "--in-fmt", dest='fmt_infile', nargs=1, type=str,
                   metavar="FORMAT", default=cli.DEFAULT_SEQ_FMT,
                   choices=cli.AVAIL_SEQ_FMTS,
                   help=("sequence file format of input"
                         " DEFAULT: {}").format(cli.DEFAULT_SEQ_FMT))
    return p

def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__,
                                     parents=[cli.get_base_parser(),
                                              cli.get_seq_out_parser(),
                                              match_order_parser(),
                                              ])
    args = parser.parse_args(argv)
    return args

def main():
    args = parse_args(sys.argv[1:])
    logging.basicConfig(level=args.log_level)
    logger.debug(args)

    for rec in order_recs(parse(args.ord_handle, args.fmt_infile),
                          index(args.match_path, args.fmt_infile)):
        write(rec, args.out_handle, args.fmt_outfile)

if __name__ == '__main__':
    main()
