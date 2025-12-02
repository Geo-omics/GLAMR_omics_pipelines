"""
Remove primer sequences from paired-end fastq-formated reads
"""
import argparse
from collections import Counter
from contextlib import ExitStack
from enum import Enum, auto
import gzip
from itertools import batched
import json
from os import environ
from pathlib import Path
import shutil
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


class Prof(Enum):
    """ primer test profile """
    STRICT = auto()
    RELAXED = auto()
    CONSISTENT = auto()


def find_5p_primer(primer, fastq, test_profile, log):
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
        data = []
        discards = []
        for idx, line in enumerate(info_file):
            row = line.rstrip('\n').split('\t')
            stats[('err', row[1])] += 1
            if row[1] == '0':
                # no error
                _, _, a, b, seq, *_ = row
                stats[(a, b)] += 1
                data.append((idx, int(a), int(b)))
            else:
                # always discard those with error
                discards.append(idx)
        if environ.get('KEEP_CUTADAPT_INFO'):
            dst = Path(fastq).with_suffix('.cutadapt_info.txt')
            shutil.copyfile(info_file.name, dst)
            print(f'[OK] Saving {dst}', file=log)

    match test_profile:
        case Prof.STRICT:
            for idx, a, b in data:
                if a == '0' and b == len(primer.sequence):
                    # keep if perfect
                    continue
                else:
                    discards.append(idx)
        case Prof.RELAXED:
            # keep all
            pass
        case Prof.CONSISTENT:
            (aa, bb), _ = Counter((a, b) for _, a, b in data).most_common()[0]
            print(f'Most common case: {aa=} {bb=}', file=log)
            for idx, a, b in data:
                if a == aa and b == bb:
                    # keep the most common case
                    continue
                else:
                    discards.append(idx)
        case _:
            raise ValueError('invalid test profile')
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

        fwd_pr, rev_pr, swapped, clean, rc = load_target_info(target_info)

        if swapped:
            raise ValueError('single data should not be swapped')

        if rc:
            # What we call the reverse primer is actually sequenced at the
            # single read's 5' end, so just for cutadapts benefit, here we're
            # swapping the notion of fwd and rev primers
            fwd_clean = clean['rev']
            rev_clean = clean['fwd']
            tmp_pr = fwd_pr
            fwd_pr = rev_pr
            rev_pr = tmp_pr
            del tmp_pr
        else:
            fwd_clean = clean['fwd']
            rev_clean = clean['rev']

        args = []  # args for cutadapt
        try:
            if not fwd_clean:
                args += ['-g', fwd_pr.sequence]
            if not rev_clean:
                args += ['-a', rev_pr.sequence]
        except AttributeError as e:
            # e.g. primer is None
            raise RuntimeError(
                f'something is marked "not clean" but no corresponding primer '
                f'was given? {e}'
            )
        if not args:
            raise RuntimeError('No args?  Something should not be clean!')

        discards = find_5p_primer(fwd_pr, single_fastq, Prof.CONSISTENT, log)
        print(f'discarding {len(discards)} read', file=log)

        with NamedTemporaryFile('w+t') as single_tmp:
            filter_fastq(discards, single_fastq, single_tmp)
            cmd = [
                'cutadapt',
                *args,
                '--discard-untrimmed',
                '--output', single_out,
                '--info-file', Path(single_out).with_name('trim_info.txt'),
                single_tmp.name,
            ]
            print('command:\n', *cmd, file=log)
            run(cmd, check=True, stdout=log)


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

        fwd_pr, rev_pr, swapped, clean, rc = load_target_info(target_info)

        if rc:
            raise ValueError('RC\'d paired data should not be a thing')

        if swapped:
            print('WARNING: swapping fwd <-> rev fastq files',
                  flush=True, file=log)
            temp = fwd_fastq
            fwd_fastq = rev_fastq
            rev_fastq = temp
            del temp

        args = []
        try:
            if not clean['fwd']:
                args += ['-g', fwd_pr.sequence]
            if not clean['fwd_rev']:
                args += ['-a', rev_pr.sequence]
            if not clean['rev']:
                args += ['-G', rev_pr.sequence]
            if not clean['rev_fwd']:
                args += ['-A', fwd_pr.sequence]
        except AttributeError as e:
            # e.g. primer is None
            raise RuntimeError(
                f'something is marked "not clean" but no corresponding '
                f'primer was given? {e}'
            )
        if not args:
            raise RuntimeError('No args?  Something should not be clean!')

        fwd_discards = find_5p_primer(fwd_pr, fwd_fastq, Prof.CONSISTENT, log)
        rev_discards = find_5p_primer(rev_pr, rev_fastq, Prof.CONSISTENT, log)
        discards = set(fwd_discards).union(rev_discards)
        print(f'discarding {len(discards)} read pairs', file=log)

        fwd_tmp = NamedTemporaryFile('w+t')
        rev_tmp = NamedTemporaryFile('w+t')
        with fwd_tmp as fwd_tmp, rev_tmp as rev_tmp:
            filter_fastq(discards, fwd_fastq, fwd_tmp)
            filter_fastq(discards, rev_fastq, rev_tmp)

            cmd = [
                'cutadapt',
                *args,
                '--discard-untrimmed',
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

    if 'swapped' in info:
        swapped = info['swapped']
    elif info['layout'] == 'single':
        swapped = False
    else:
        raise ValueError('info does not have "swapped" key')

    if not isinstance(swapped, bool):
        raise ValueError(f'swapped must be boolean: {swapped}')

    rc = info.get('model_rc', False)

    clean = dict(
        fwd=info.get('fwd_clean'),
        rev=info.get('rev_clean'),
        fwd_rev=info.get('fwd_rev_clean'),
        rev_fwd=info.get('rev_fwd_clean'),
    )
    return fwd_primer, rev_primer, swapped, clean, rc


if __name__ == '__main__':
    cli()
