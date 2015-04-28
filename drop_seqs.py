#!/usr/bin/env python3
"""Given a list of sequence IDs and a sequence file, output the latter
without any of the former.

"""

from Bio.SeqIO import parse, write
import sys
import argparse
import cli
import logging

logger = logging.getLogger(__name__)

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

def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__,
                                     parents=[cli.get_base_parser(),
                                              cli.get_seq_out_parser(),
                                              cli.get_list_in_parser(),
                                              cli.get_seq_in_parser(),
                                              ])
    args = parser.parse_args(argv[1:])
    return args

def main():
    args = parse_args(sys.argv)
    logging.basicConfig(level=args.log_level)
    logger.debug(args)

    for rec in rm_recs(parse(args.in_handle, args.fmt_infile),
                       get_list(args.list_handle)):
        write(rec, args.out_handle, args.fmt_outfile)

if __name__ == '__main__':
    main()
