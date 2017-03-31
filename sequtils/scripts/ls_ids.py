#!/usr/bin/env python3
"""Pull IDs from sequence files.

"""

from Bio.SeqIO import parse
import sys
import argparse
from ..lib import cli
from signal import signal, SIGPIPE, SIG_DFL
import logging

logger = logging.getLogger(__name__)
logging.captureWarnings(True)

def get_ids(recs):
    for rec in recs:
        yield rec.id

def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__,
                                     parents=[cli.get_base_parser(),
                                              cli.get_seq_in_parser(),
                                              ])
    args = parser.parse_args(argv[1:])
    return args

def main():
    args = parse_args(sys.argv)
    logging.basicConfig(level=args.log_level)
    logger.debug(args)

    id_gen = get_ids(parse(args.in_handle, args.fmt_infile))
    for id in id_gen:
        args.out_handle.write("{}\n".format(id))

if __name__ == "__main__":
    signal(SIGPIPE,SIG_DFL)
    main()
