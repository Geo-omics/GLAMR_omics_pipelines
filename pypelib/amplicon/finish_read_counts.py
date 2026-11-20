#!/usr/bin/env python3
"""
Finalize the read count stats table from dada2

Outputs to stdout.
"""
import argparse
from pathlib import Path
from sys import stderr


argp = argparse.ArgumentParser(description=__doc__)
argp.add_argument(
    'track_read_counts',
    help='The track_read_counts.tsv file as made by the dada2 script',
)
argp.add_argument(
    'raw_reads_stats',
    nargs='*',
    help='The raw reads stats, one per sample, as made by seqkit',
)
args = argp.parse_args()


def msg(*args, **kwargs):
    print(*args, **kwargs, file=stderr)


msg(f'Reading {args.track_read_counts} ...', end='', flush=True)
with open(args.track_read_counts) as ifile:
    header = ifile.readline().split()
    rows = [line.split() for line in ifile]
msg(f'{len(rows)} rows [OK]')

raw_stats_files0 = set(args.raw_reads_stats)
if len(raw_stats_files0) != len(rows):
    msg(f'[WARNING]: {args.track_read_counts} has data for {len(rows)} samples'
        f'but {len(raw_stats_files0)} distinct raw read stat files are '
        f'provided')


# Check raw stats paths and see what sample each belongs to
raw_stats_files = {}
for row in rows:
    sample_id = row[0]
    for path in raw_stats_files0:
        if sample_id in Path(path).parts:
            raw_stats_files[sample_id] = path
            break
    else:
        argp.error(f'No raw reads stats for sample {sample_id} was provided')


msg('Reading raw reads stats... ', end='', flush=True)
raw_data = {}
for sample_id, path in raw_stats_files.items():
    with open(path) as ifile:
        try:
            num_seqs_idx = ifile.readline().split().index('num_seqs')
        except ValueError:
            argp.error(f'failed parsing header of {path}')

        try:
            nums = [
                int(line.strip().split('\t')[num_seqs_idx])
                for line in ifile
            ]
        except (IndexError, ValueError) as e:
            argp.error(f'Failed parsing {path}: {e}')

        match len(nums):
            case 1:
                # single-ended data
                raw_data[sample_id] = nums.pop()
            case 2:
                # paired-end data
                a, b = nums
                if a != b:
                    argp.error(f'fwd-vs-rev raw read count differs in {path}: {nums=}')

                raw_data[sample_id] = a
            case _:
                argp.error(f'expecting two rows: {path}')
msg('[OK]')

# Output
print(header[0], 'raw', *header[1:], sep='\t')
for row in rows:
    sample_id = row[0]
    print(sample_id, raw_data[sample_id], *row[1:], sep='\t')
