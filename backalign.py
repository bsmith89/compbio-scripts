#!/usr/bin/env python3
"""Align codons based on a corresponding amino-acid sequence.

Not defensive in the slightest:
    Doesn't check that codons and AAs match
    Doesn't check that sequence labels match (see "match_order.py")

"""

from Bio.SeqIO import parse, write
from Bio.Seq import Seq
import sys
import argparse
from copy import copy

DEFAULT_INFILE_FMT = 'fasta'
DEFAULT_OUTFILE_FMT = 'fasta'

def codons(sequence):
    """Return groups of 3 items from a sequence.
    
    This is specifically meant to return codons from an in-frame sequence."""
    codon = ""
    for nucl in sequence:
        codon += nucl
        if len(codon) == 3:
            yield codon
            codon = ""

def backalign(nucl, prot):
    align_nucl = ""
    codon_iter = codons(nucl)
    for amino in prot:
        if amino == ".":
            continue
        elif amino.islower():
            next(codon_iter)
        elif amino == "-":
            align_nucl += "---"
        else:
            align_nucl += next(codon_iter)
    return Seq(align_nucl)

def backalign_recs(nucl_recs, prot_recs):
    for nucl_rec, prot_rec in zip(nucl_recs, prot_recs):
        out_rec = copy(nucl_rec)
        out_rec.seq = backalign(nucl_rec.seq, prot_rec.seq)
        yield out_rec

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
    p.add_argument('in_prot', type=argparse.FileType('r'),
                   metavar="PROT-FILE")
    p.add_argument('in_nucl', nargs='?', type=argparse.FileType('r'),
                   metavar="NUCL-FILE",
                   default=sys.stdin)

    args = p.parse_args()

    for rec in backalign_recs(parse(args.in_nucl, args.fmt_infile),
                              parse(args.in_prot, args.fmt_infile)):
        write(rec, sys.stdout, args.fmt_outfile)

if __name__ == '__main__':
    main()
