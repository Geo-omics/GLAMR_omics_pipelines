#!/usr/bin/env python3
"""
Fix the taxonomy output file by substituting the first column sequences with IDs.
"""
import argparse
from contextlib import ExitStack
from sys import stdout

from Bio import SeqIO


def cli():
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument(
        'asv_fasta',
        help='Fasta file with ASV sequences as made by the amplicon_asv_check rule',
    )
    argp.add_argument(
        'dada2_tax_assignment',
        help='taxonomy assignment file made with the assign_taxonomy.R script.'
    )
    argp.add_argument(
        '--output', '-o',
        help='Output file.  Print to stdout if omitted.'
    )
    args = argp.parse_args()
    main(args.asv_fasta, args.dada2_tax_assignment, args.output)


def main(fasta, taxonomy, output):
    seq2id = {str(i.seq): i.id for i in SeqIO.parse(fasta, 'fasta')}
    lines = []
    with open(taxonomy) as ifile:
        lines.append(ifile.readline())  # header
        for lnum, line in enumerate(ifile, start=2):
            seq, _, rest = line.partition('\t')
            try:
                asv_id = seq2id[seq]
            except KeyError:
                raise LookupError(f'sequence on line {lnum} not in fasta file')
            lines.append(asv_id + '\t' + rest)

    with ExitStack() as estack:
        if output:
            ofile = estack.enter_context(open(output, 'w'))
        else:
            ofile = stdout
        ofile.writelines(lines)


if __name__ == '__main__':
    cli()
