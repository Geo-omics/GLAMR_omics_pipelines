"""
Remove primer sequences from paired-end fastq-formated reads
"""
import argparse
from collections import Counter
from contextlib import ExitStack
import gzip
from itertools import batched
import json
from pathlib import Path
from subprocess import run, DEVNULL
from tempfile import NamedTemporaryFile

from . import get_models


def cli():
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument('target_info', help='The target info file.')
    argp.add_argument(
        'raw_fwd_fastq',
        help='Forward (or single) raw fastq file, gzipped',
    )
    argp.add_argument(
        'raw_rev_fastq',
        required=False,
        help='Reverse raw fastq file, gzipped',
    )
    argp.add_argument('--fwd-out', required=True,
                      help='forward (or single) output file')
    argp.add_argument('--rev-out', help='reverse output file')
    argp.add_argument('--log', help='log file, will otherwise go to stdout')
    args = argp.parse_args()

    if args.raw_rev_fastq:
        main_paired(
            target_info=args.target_info,
            fwd_fastq=args.raw_fwd_fastq,
            rev_fastq=args.raw_rev_fastq,
            fwd_out=args.fwd_out,
            rev_out=args.rev_out,
            log=args.log,
        )
    else:
        main_single(
            target_info=args.target_info,
            single_fastq=args.raw_fwd_fastq,
            single_out=args.fwd_out,
            log=args.log,
        )


def main_single(
    target_info=None,
    single_fastq=None,
    single_out=None,
    log=None,
):
    with ExitStack() as estack:
        if log:
            log = estack.enter_context(open(log, 'w'))
        else:
            log = None

        fwd_primer, rev_primer, swapped, clean = load_target_info(target_info)

        if swapped:
            raise ValueError('single data should not be swapped')

        args = []
        try:
            if not clean['fwd']:
                args += ['-g', fwd_primer.sequence]
            if not clean['rev']:
                args += ['-a', rev_primer.sequence]
        except AttributeError as e:
            # e.g. primer is None
            raise RuntimeError(
                f'something is marked "not clean" but no corresponding primer '
                f'was given? {e}'
            )
        if not args:
            raise RuntimeError('No args?  Something should not be clean!')

        cmd = [
            'cutadapt',
            *args,
            '--discard-untrimmed',
            '--output', single_out,
            '--info-file', Path(single_out).with_name('trim_info.txt'),
            single_fastq,
        ]
        print('command:\n', *cmd, file=log)
        run(cmd, check=True, stdout=log)


def find_5p_primer(primer, fastq, log):
    """
    Get indexes of sequences with imperfect 5' primer
    """
    stats = Counter()
    with NamedTemporaryFile('w+t') as info_file:
        cmd = [
            'cutadapt',
            '-g', primer.sequence,
            '--info-file', info_file.name,
            fastq,
        ]
        run(cmd, check=True, stdout=DEVNULL, stderr=log)
        info_file.seek(0)
        discards = []
        for idx, line in enumerate(info_file):
            row = line.rstrip('\n').split('\t')
            stats[('err', row[1])] += 1
            if row[1] == '0':
                # no error
                _, _, a, b, seq, *_ = row
                stats[(a, b)] += 1
                if a == '0' and int(b) == len(primer.sequence):
                    # keep
                    continue
            discards.append(idx)
    good_count = idx + 1 - len(discards)
    print(
        f'{fastq}: keeping {good_count} reads, discarding '
        f'{len(discards) / (idx + 1):.1%}\nstats: {stats}',
        file=log
    )
    return discards


def filter_fastq(discards, infilename, outfile):
    with gzip.open(infilename, 'rt') as ifile:
        for idx, (head, seq, plus, qual) in enumerate(batched(ifile, n=4)):
            if not head.startswith('@'):
                raise ValueError('not a fastq header')
            if not plus.startswith('+'):
                raise ValueError('not a plus line')
            if idx in discards:
                continue
            outfile.write(head)
            outfile.write(seq)
            outfile.write(plus)
            outfile.write(qual)
    outfile.flush()


def main_paired(
    target_info=None,
    fwd_fastq=None,
    rev_fastq=None,
    fwd_out=None,
    rev_out=None,
    log=None,
):
    with ExitStack() as estack:
        if log:
            log = estack.enter_context(open(log, 'w'))
        else:
            log = None

        fwd_primer, rev_primer, swapped, clean = load_target_info(target_info)

        if swapped:
            print('WARNING: swapping fwd <-> rev fastq files',
                  flush=True, file=log)
            temp = fwd_fastq
            fwd_fastq = rev_fastq
            rev_fastq = temp
            del temp

        fwd_discards = find_5p_primer(fwd_primer, fwd_fastq, log)
        rev_discards = find_5p_primer(rev_primer, rev_fastq, log)
        discards = set(fwd_discards).union(rev_discards)
        print(f'discarding {len(discards)} read pairs', file=log)

        fwd_tmp = NamedTemporaryFile('w+t')
        rev_tmp = NamedTemporaryFile('w+t')
        with fwd_tmp as fwd_tmp, rev_tmp as rev_tmp:
            filter_fastq(discards, fwd_fastq, fwd_tmp)
            filter_fastq(discards, rev_fastq, rev_tmp)

            args = []
            try:
                if not clean['fwd']:
                    args += ['-g', fwd_primer.sequence]
                if not clean['fwd_rev']:
                    args += ['-a', rev_primer.sequence]
                if not clean['rev']:
                    args += ['-G', rev_primer.sequence]
                if not clean['rev_fwd']:
                    args += ['-A', fwd_primer.sequence]
            except AttributeError as e:
                # e.g. primer is None
                raise RuntimeError(
                    f'something is marked "not clean" but no corresponding '
                    f'primer was given? {e}'
                )
            if not args:
                raise RuntimeError('No args?  Something should not be clean!')

            cmd = [
                'cutadapt',
                *args,
                '--discard-untrimmed',
                # '-u', '17', '-U', '20',
                '--output', fwd_out,
                '--paired-output', rev_out,
                '--info-file', Path(fwd_out).with_name('trim_info.txt'),
                fwd_tmp.name,
                rev_tmp.name,
            ]
            print('command:\n', *cmd, file=log)
            run(cmd, check=True, stdout=log)


def main_paired_gen_1(
    target_info=None,
    fwd_fastq=None,
    rev_fastq=None,
    fwd_out=None,
    rev_out=None,
    log=None,
):
    """ DEPRECATED """
    with ExitStack() as estack:
        if log:
            log = estack.enter_context(open(log, 'w'))
        else:
            log = None

        fwd_primer, rev_primer, swapped, clean = load_target_info(target_info)

        if swapped:
            print('WARNING: swapping fwd <-> rev fastq files',
                  flush=True, file=log)

        args = []
        try:
            if not clean['fwd']:
                args += ['-g', fwd_primer.sequence]
            if not clean['fwd_rev']:
                args += ['-a', rev_primer.sequence]
            if not clean['rev']:
                args += ['-G', rev_primer.sequence]
            if not clean['rev_fwd']:
                args += ['-A', fwd_primer.sequence]
        except AttributeError as e:
            # e.g. primer is None
            raise RuntimeError(
                f'something is marked "not clean" but no corresponding primer '
                f'was given? {e}'
            )
        if not args:
            raise RuntimeError('No args?  Something should not be clean!')

        cmd = [
            'cutadapt',
            *args,
            '--discard-untrimmed',
            # '-u', '17', '-U', '20',
            '--output', fwd_out,
            '--paired-output', rev_out,
            '--info-file', Path(fwd_out).with_name('trim_info.txt'),
            rev_fastq if swapped else fwd_fastq,
            fwd_fastq if swapped else rev_fastq,
        ]
        print('command:\n', *cmd, file=log)
        run(cmd, check=True, stdout=log)


def load_target_info(target_info):
    with open(target_info) as ifile:
        info = json.load(ifile)

    if not (model_name := info.get('model_name')):
        raise ValueError('model_name is required')

    try:
        model = get_models()[model_name]
    except KeyError as e:
        raise ValueError(f'unknown model: {model_name}') from e

    if fwd_primer := info.get('fwd_primer'):
        for i in model.fwd_primers:
            if i.name == fwd_primer:
                fwd_primer = i
                break
        else:
            raise ValueError(
                f'Not a fwd primer in model {model_name}: {fwd_primer}'
            )

    if rev_primer := info.get('rev_primer'):
        for i in model.rev_primers:
            if i.name == rev_primer:
                rev_primer = i
                break
        else:
            raise ValueError(
                f'Not a rev primer in model {model_name}: {rev_primer}')

    swapped = info['swapped']
    if not isinstance(swapped, bool):
        raise ValueError(f'swapped must be boolean: {swapped}')

    clean = dict(
        fwd=info.get('fwd_clean'),
        rev=info.get('rev_clean'),
        fwd_rev=info.get('fwd_rev_clean'),
        rev_fwd=info.get('rev_fwd_clean'),
    )
    return fwd_primer, rev_primer, swapped, clean


if __name__ == '__main__':
    cli()
