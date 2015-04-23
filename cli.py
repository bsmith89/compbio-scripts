"""Common command-line interface elements for sequence analysis utilities.

"""

import argparse
import Bio

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
                   help=("file format of infile or stdin"
                         " DEFAULT: {}").format(DEFAULT_SEQ_FMT))
    p.add_argument("-t", "--out-fmt", dest='fmt_outfile', nargs=1, type=str,
                   metavar="FORMAT", default=DEFAULT_SEQ_FMT,
                   choices=AVAIL_SEQ_FMTS,
                   help=("file format of outfile or stdout"
                         " DEFAULT: {}").format(DEFAULT_SEQ_FMT))
    p.add_argument("-v", "--verbose",
                   dest='log_level', action='store_const', const=10,
                   help=("set loggin level to 10 (debug)"))
    p.add_argument("--log-level", type=int,
                   default=DEFAULT_LOG_LVL,
                   help=("logging level (higher=fewer messages)"
                         " DEFAULT: {}").format(DEFAULT_LOG_LVL))
    return p

if __name__ == "__main__":
    pass
