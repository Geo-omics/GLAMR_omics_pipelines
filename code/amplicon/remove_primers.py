"""
Remove primer sequences from paired-end fastq-formated reads
"""
import argparse
from subprocess import run

from .primers import Primer


def main():
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument('target_info', help='The target info file.')
    argp.add_argument('raw_fwd_fastq', help='Forward raw fastq file, gzipped')
    argp.add_argument('raw_rev_fastq', help='Reverse raw fastq file, gzipped')
    argp.add_argument('--fwd-out', help='forward output file')
    argp.add_argument('--rev-out', help='reverse output file')
    args = argp.parse_args()

    fprim, rprim, swapped = get_primer_seqs(args.target_info)
    run_cutadapt(fprim, rprim, swapped, args)


def get_primer_seqs(target_info):
    with open(target_info) as ifile:
        f = r = None
        swapped = None
        for line in ifile:
            line = line.strip()
            if line.startswith('primers:'):
                value = line.removeprefix('primers:').strip()
                f, _, r = value.partition('/')
            elif line.startswith('swapped:'):
                swapped = line.removeprefix('swapped:').strip()

    if not f and not r:
        raise RuntimeError('no primer information found')

    fwd_primer = rev_primer = None
    for i in Primer.load():
        if i.name == f:
            if i.direction != Primer.FWD:
                raise RuntimeError(f'wrong direction: {f} vs. {i}')
            fwd_primer = i
        elif i.name == r:
            if i.direction != Primer.REV:
                raise RuntimeError(f'wrong direction: {r} vs. {i}')
            rev_primer = i
    if fwd_primer is None:
        raise RuntimeError(f'fwd primer not found: {f}')
    if rev_primer is None:
        raise RuntimeError(f'rev primer not found: {r}')

    if swapped is None:
        swapped = False
    elif swapped == 'YES':
        swapped = True
    else:
        raise ValueError(f'invalid swapped value: {swapped}')

    return fwd_primer, rev_primer, swapped


def run_cutadapt(fwd_primer, rev_primer, swapped, args):
    if swapped:
        print('WARNING: swapping fwd <-> rev fastq files', flush=True)
    cmd = [
        'cutadapt',
        '-g', fwd_primer.sequence,
        '-G', rev_primer.sequence,
        '-o', args.fwd_out,
        '-p', args.rev_out,
        args.raw_rev_fastq if swapped else args.raw_fwd_fastq,
        args.raw_fwd_fastq if swapped else args.raw_rev_fastq,
    ]
    run(cmd, check=True)


if __name__ == '__main__':
    main()
