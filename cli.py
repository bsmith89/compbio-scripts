"""Common command-line interface elements for sequence analysis utilities.

"""

import argparse
import Bio
import sys

DEFAULT_SEQ_FMT = 'fasta'
# TODO: Provide more options here.  Presumably anything available in Biopython.
AVAIL_SEQ_FMTS = ['fasta', 'fastq', 'genbank', 'sff', 'swiss', 'tab']
DEFAULT_LOG_LVL = 30

def get_default_parser():
    """Get a command-line argument parser with standard elements.

    This ArgumentParser is meant to be passed as a parent to a local
    ArgumentParser initialization.

    """
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("-f", "--in-fmt", dest='fmt_infile', nargs=1, type=str,
                   metavar="FORMAT", default=DEFAULT_SEQ_FMT,
                   choices=AVAIL_SEQ_FMTS,
                   help=("file format of input"
                         " DEFAULT: {}").format(DEFAULT_SEQ_FMT))
    p.add_argument("-t", "--out-fmt", dest='fmt_outfile', nargs=1, type=str,
                   metavar="FORMAT", default=DEFAULT_SEQ_FMT,
                   choices=AVAIL_SEQ_FMTS,
                   help=("file format of output"
                         " DEFAULT: {}").format(DEFAULT_SEQ_FMT))
    p.add_argument("--log-level", type=int,
                   default=DEFAULT_LOG_LVL,
                   help=("logging level (higher=fewer messages)"
                         " DEFAULT: {}").format(DEFAULT_LOG_LVL))
    p.add_argument("-v", "--verbose",
                   dest='log_level', action='store_const', const=10,
                   help=("set loggin level to 10 (debug)"))
    p.add_argument("-o", "--out-file", dest='out_handle',
                   type=argparse.FileType('w'),
                   metavar="OUTFILE", default=sys.stdout,
                   help=("write to file instead of stdout"))
    return p

def get_infile_parser():
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("in_handle", nargs='?', type=argparse.FileType('r'),
                   metavar="INFILE", default=sys.stdin,
                   help=("sequences"))
    return p

if __name__ == "__main__":
    pass
