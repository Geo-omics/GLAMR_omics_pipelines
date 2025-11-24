"""
Set up test environment directory structure
"""
import argparse
from collections import Counter
from pathlib import Path
import shutil

from .utils import UsageError


def cli():
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument(
        'source_directory',
        help='Path to existing (production?) data storage, usually the GLAMR '
             'data/ directory',
    )
    argp.add_argument(
        'destination',
        help='Destination directory (the parent of the test environment, so '
             'without the usual "data".)',
    )
    argp.add_argument(
        '--sample-type-only',
        help='Only process given sample type (usually the directory name under'
             'data/omics)'
    )
    argp.add_argument(
        '--exist-ok',
        action='store_true',
        help='Do not raise an error on exsting files and directories',
    )
    argp.add_argument(
        '--copy-raw-reads',
        action='store_true',
        help='Optionally copy any raw reads file too',
    )
    argp.add_argument(
        '--traceback',
        action='store_true',
        help='print full traceback on all errors',
    )
    args = argp.parse_args()
    try:
        main(
            args.source_directory,
            args.destination,
            exist_ok=args.exist_ok,
            raw_reads=args.copy_raw_reads,
            sample_type=args.sample_type_only,
        )
    except UsageError as e:
        if args.traceback:
            raise
        argp.error(e)


def main(src, dst, exist_ok=False, raw_reads=False, sample_type=None):
    src = Path(src)
    if not src.is_dir():
        raise UsageError(f'no such directory: {src}')

    dst = Path(dst)
    if not src.is_dir():
        raise UsageError(f'no such directory: {dst}')

    data = Path(src.name)
    src = src.parent  # use as src root from now on
    omics = data / 'omics'
    projects = data / 'projects'

    sample_types = [
        i.name
        for i in (src / omics).iterdir()
        if i.is_dir() and sample_type in (None, i.name)
    ]
    if not sample_types:
        raise UsageError(
            'invalid sample type passed or there are no samp type dirs at all'
        )

    (dst / data).mkdir(exist_ok=exist_ok)
    (dst / omics).mkdir(exist_ok=exist_ok)
    stats = Counter()
    for stype in sample_types:
        stype_dir = omics / stype
        (dst / stype_dir).mkdir(exist_ok=exist_ok)
        for i in (src / stype_dir).glob('samp_*'):
            if not i.is_dir():
                continue
            sample_dir = stype_dir / i.name
            (dst / sample_dir).mkdir(exist_ok=exist_ok)
            stats[(stype,)] += 1
            reads = sample_dir / 'reads'
            accn = src / reads / 'accession'
            if accn.is_file():
                (dst / reads).mkdir(exist_ok=exist_ok)
                shutil.copy(accn, dst / reads)
                stats[(stype, 'accn')] += 1
            if raw_reads:
                raws = [
                    src / reads / 'raw_single_reads.fastq.gz',
                    src / reads / 'raw_fwd_reads.fastq.gz',
                    src / reads / 'raw_rev_reads.fastq.gz',
                ]
                for fq in raws:
                    if fq.is_file():
                        shutil.copy(fq, dst / reads)
                        stats[(stype, 'raw_' + fq.name.split('_')[1])] += 1

    (dst / projects).mkdir(exist_ok=exist_ok)
    for i in (src / projects).glob('set_*'):
        if not i.is_dir():
            continue
        dset_dir = projects / i.name
        (dst / dset_dir).mkdir(exist_ok=exist_ok)
        for stype in sample_types:
            if not (src / dset_dir / stype).is_dir():
                continue
            (dst / dset_dir / stype).mkdir(exist_ok=exist_ok)
            for j in (src / dset_dir / stype).glob('samp_*'):
                if j.is_symlink():
                    try:
                        (dst / dset_dir / stype / j.name).symlink_to(
                            dst / omics / stype / j.name
                        )
                    except FileExistsError:
                        if not exist_ok:
                            raise
                    stats[(i.name, stype)] += 1

    for key, count in sorted(stats.items()):
        key = ' // '.join(str(i) for i in key)
        print(f'  {key:<30}: {count:>5}')


if __name__ == '__main__':
    cli()
