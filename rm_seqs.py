#!/usr/bin/env python3
"""Given a list of sequence IDs and a sequence file, output the latter
without any of the former.

"""

from Bio.SeqIO import parse, write
import sys
import argparse

DEFAULT_INFILE_FMT = 'fasta'
DEFAULT_OUTFILE_FMT = 'fasta'

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
    p = argparse.ArgumentParser(description=__doc__)
    # Arguments which may be generalizable to many scripts.
    # TODO: Break these out into a module which can be imported.
    p.add_argument("--in-fmt", "-f", dest='fmt_infile', nargs=1, type=str,
                   metavar="FORMAT", default=DEFAULT_INFILE_FMT,
                   help="file format of infile or stdin")
    p.add_argument("--out-fmt", "-t", dest='fmt_outfile', nargs=1, type=str,
                   metavar="FORMAT", default=DEFAULT_OUTFILE_FMT,
                   help="file format of outfile or stdout")
    p.add_argument("--verbose", "-v", action='count',
                   help="increase verbosity")

    # Arguments specific to this script.
    p.add_argument('rm_handle', type=argparse.FileType('r'),
                   metavar="LISTFILE")
    p.add_argument('in_handle', nargs='?', type=argparse.FileType('r'),
                   metavar="SEQFILE", default=sys.stdin)

    args = p.parse_args()

    for rec in rm_recs(parse(args.in_handle, args.fmt_infile),
                       get_list(args.rm_handle)):
        write(rec, sys.stdout, args.fmt_outfile)

if __name__ == '__main__':
    main()
