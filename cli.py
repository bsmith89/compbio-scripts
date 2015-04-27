"""Common command-line interface elements for sequence analysis utilities.

"""

import argparse
import sys
import logging


INPUT_GROUP = "input arguments", "format and origin"
OUTPUT_GROUP = "output arguments", "format and destination"
LOG_GROUP = "logging arguments", "verbosity and debugging"

# TODO: Provide more options here.  Presumably anything available in Biopython.
AVAIL_SEQ_FMTS = ['fasta', 'fastq', 'genbank', 'sff', 'swiss', 'tab']
DEFAULT_SEQ_FMT = 'fasta'

AVAIL_ALIGN_FMTS = ['clustal', 'emboss', 'fasta', 'nexus', 'phylip',
                    'phylip-sequential', 'phylip-relaxed', 'stockholm']
DEFAULT_ALIGN_FMT = 'fasta'

DEFAULT_LOG_LVL = 30

def get_base_parser():
    p = argparse.ArgumentParser(add_help=False)
    g = p.add_argument_group(*LOG_GROUP)
    g.add_argument("--log-level", type=int,
                   default=DEFAULT_LOG_LVL,
                   help=("logging level (higher=fewer messages)"
                         " DEFAULT: {}").format(DEFAULT_LOG_LVL))
    g.add_argument("-v", "--verbose",
                   dest='log_level', action='store_const', const=10,
                   help=("set loggin level to 10 (debug)"))
    g.add_argument("-V", "--version", dest='show_version',
                   action='store_true',
                   help=("show version information and exit"))
    h = p.add_argument_group(*OUTPUT_GROUP)
    h.add_argument("-o", "--out-file", dest='out_handle',
                   type=argparse.FileType('w'),
                   metavar="OUTFILE", default=sys.stdout,
                   help=("write to file instead of stdout"))
    return p

def get_seq_in_parser(optional=True):
    if optional:
        in_handle_nargs = '?'
    else:
        in_handle_nargs = None

    p = argparse.ArgumentParser(add_help=False)
    g = p.add_argument_group(*INPUT_GROUP)
    g.add_argument("in_handle", nargs=in_handle_nargs,
                   type=argparse.FileType('r'),
                   metavar="SEQUENCE", default=sys.stdin,
                   help=("sequence file"))
    g.add_argument("-f", "--in-fmt", dest='fmt_infile', nargs=1, type=str,
                   metavar="FORMAT", default=DEFAULT_SEQ_FMT,
                   choices=AVAIL_SEQ_FMTS,
                   help=("sequence file format of input"
                         " DEFAULT: {}").format(DEFAULT_SEQ_FMT))
    return p


def get_seq_out_parser():
    p = argparse.ArgumentParser(add_help=False)
    g = p.add_argument_group(*OUTPUT_GROUP)
    g.add_argument("-t", "--out-fmt", dest='fmt_outfile', nargs=1, type=str,
                   metavar="FORMAT", default=DEFAULT_SEQ_FMT,
                   choices=AVAIL_SEQ_FMTS,
                   help=("sequence file format of output"
                         " DEFAULT: {}").format(DEFAULT_SEQ_FMT))
    return p

def get_align_in_parser(optional=False):
    if optional:
        in_handle_nargs = '?'
    else:
        in_handle_nargs = None

    p = argparse.ArgumentParser(add_help=False)
    g = p.add_argument_group(*INPUT_GROUP)
    g.add_argument('align_handle', nargs=in_handle_nargs,
                   type=argparse.FileType('r'),
                   metavar="ALIGNMENT",
                   help=("alignment file"))
    g.add_argument('--align-fmt', dest='fmt_align', nargs=1, type=str,
                   metavar="FORMAT", default=DEFAULT_ALIGN_FMT,
                   choices=AVAIL_ALIGN_FMTS,
                   help=("file format of aligned protein sequences"
                         " DEFAULT: {}").format(DEFAULT_ALIGN_FMT))
    return p

def get_list_in_parser():
    p = argparse.ArgumentParser(add_help=False)
    g = p.add_argument_group(*INPUT_GROUP)
    g.add_argument('list_handle', help="list of sequence IDs",
                   metavar='LIST', type=argparse.FileType('r'))
    return p

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Test parser.",
                                     parents=[get_base_parser(),
                                              get_seq_in_parser(),
                                              get_align_in_parser(),
                                              get_seq_out_parser()])
    args = parser.parse_args(argv)

    logger = logging.getLogger(__name__)
    logging.basicConfig(level=args.log_level)
    logger.debug(args)

    return args

def main():
    args = parse_args(sys.argv[1:])

if __name__ == "__main__":
    main()
