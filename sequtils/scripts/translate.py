#!/usr/bin/env python3
"""Translate nucleotide sequences to amino acids."""

from Bio.SeqIO import parse, write
import sys
import argparse
from copy import copy
from ..lib import cli
import logging
from signal import signal, SIGPIPE, SIG_DFL

logger = logging.getLogger(__name__)
logging.captureWarnings(True)

DEFAULT_CODE = 1

def translate_recs(records, code):
    """Translate sequence from records.

    This geneator takes a iterable of Bio.SeqRecord objects and the same
    records with sequences translated from nucleotides to amino acids.

    """
    for rec in records:
        # TODO: Make sure that this shallow copy of rec
        # is safe.  Is rec.seq a mutable sequence?
        # No, right?
        rec = copy(rec)
        rec.seq = rec.seq.translate(table=code)
        yield rec

def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__,
                                     parents=[cli.get_base_parser(),
                                              cli.get_seq_out_parser(),
                                              cli.get_seq_in_parser(),
                                              ])
    g = parser.add_argument_group(*cli.PAR_GROUP)
    g.add_argument('--code', type=int,
                   help=("NCBI translation table number "
                         "DEFAULT: {}").format(DEFAULT_CODE),
                   default=DEFAULT_CODE)
    args = parser.parse_args(argv[1:])
    return args

def main():
    signal(SIGPIPE,SIG_DFL)
    args = parse_args(sys.argv)
    logging.basicConfig(level=args.log_level)
    logger.debug(args)

    for trans_rec in translate_recs(parse(args.in_handle, args.fmt_infile), code=args.code):
        write(trans_rec, args.out_handle, args.fmt_outfile)

if __name__ == '__main__':
    main()
