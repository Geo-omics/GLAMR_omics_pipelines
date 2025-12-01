#!/usr/bin/env python3
"""
Guess the amplicon target based on summaries of hmmscan

Output is written to stdout.
"""
import argparse
import json
from pathlib import Path

from ..utils import load_stats, UsageError
from . import get_models, Err


def cli():
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument(
        'summary',
        nargs='+',
        help='One or more paths to summarized nhmmscan results, json formatted'
             ', as made by the amplicon.hmm_summarize script.  For paired-end '
             'sequencing two files should be provided, the first being for the'
             ' forward reads and the second for the reverse reads.',
    )
    argp.add_argument(
        '--stats',
        help='fastq statistics file as made with "seqkit stats -a" from raw '
             'fastq files',
    )
    argp.add_argument(
        '--output',
        help='Name of output file.  By default write to stdout',
    )
    argp.add_argument(
        '--traceback',
        action='store_true',
        help='always show traceback on errors',
    )

    args = argp.parse_args()
    try:
        main(args.summary, args.stats, args.output)
    except UsageError as e:
        if args.traceback:
            raise
        else:
            argp.error(e)


def main(summaries, stats_file, outfile=None):
    if isinstance(summaries, str):
        # e.g. invoked via snakemake, single file
        summaries = [summaries]

    if len(summaries) > 2:
        raise UsageError(
            f'too many items in list argument, max two summary files: '
            f'{summaries=}'
        )

    try:
        data = [read_summary(i) for i in summaries]
    except Exception as e:
        raise UsageError(
            f'failed reading summary file: {e.__class__.__name__}: {e}'
        ) from e

    try:
        lengths = get_lengths(stats_file, *summaries)
    except BadData as e:
        msg, errno = e.args
        raise UsageError(
            f'failed getting read statistic:\n'
            f'E{errno}: {msg}'
        )

    if len(data) != len(lengths):
        raise UsageError(
            'number of summary files and rows in stats file do not agree'
        )

    for dat, length in zip(data, lengths, strict=True):
        dat['length'] = length
    del lengths

    out_data = {}
    errors = []

    if len(data) == 1:
        try:
            errors = basic_checks(*data, 'single')
        except (BadData, InsufficientData) as e:
            errors.append(e.arg)
        else:
            out_data.update(check_single(*data))
    elif len(data) == 2:
        try:
            errors = basic_checks(data[0], 'fwd')
            errors += basic_checks(data[1], 'rev')
        except (BadData, InsufficientData) as e:
            errors.append(e.arg)
        else:
            out_data.update(check_paired(*data))
    else:
        raise RuntimeError('bug')

    out_data['errors'] += [
        f'E{err}: {msg}'
        for msg, err in errors
    ]
    if not out_data['errors']:
        del out_data['errors']

    out_txt = json.dumps(out_data, indent=4)
    if outfile is None:
        print(out_txt)
    else:
        with open(outfile, 'w') as ofile:
            ofile.write(out_txt)
            ofile.write('\n')
        print(f'Output written to: {outfile}')


class BadData(Exception):
    pass


class InsufficientData(Exception):
    pass


def error(errno, msg):
    print(f'error E{errno.value}:', msg)


def read_summary(file):
    with open(file) as ifile:
        data = json.load(ifile)

    all_models = get_models()
    if model_name := data.get('model'):
        model = all_models[model_name]
        data['model'] = model
    for direc, primers in [('fwd', model.fwd_primers), ('rev', model.rev_primers)]:  # noqa: E501
        key = f'{direc}_primer'
        if primer_name := data.get(key):
            for i in primers:
                if i.name == primer_name:
                    data[key] = i
                    break
            else:
                raise ValueError(
                    f'no {direc} primer with name {primer_name} in model '
                    f'{model_name}'
                )
    return data


def get_lengths(stats_file, *summaries):
    """
    Get length(s) from stats file

    Assume first row is forward data and second row i(if any) reverse.  Get the
    median (truncated towards zero) as nominal read length.
    """
    try:
        stats = load_stats(stats_file)
    except Exception as e:
        raise BadData(
            f'failed reading stats file: {e.__class__.__name__}: {e}',
            Err.E15
        ) from e

    # sanity checks
    if len(summaries) == 1 and len(stats) == 1:
        # single-ended
        pass
    elif len(summaries) == 2 and len(stats) == 2:
        # paired-end, check fwd is fwd and rev is rev
        raw_name_1 = Path(list(stats)[0]).name
        raw_name_2 = Path(list(stats)[1]).name
        sum_name_1 = Path(summaries[0]).name
        sum_name_2 = Path(summaries[1]).name
        test_raw = 'fwd' in raw_name_1 and 'rev' in raw_name_2
        test_sum = 'fwd' in sum_name_1 and 'rev' in sum_name_2
        if test_raw and test_sum:
            pass
        else:
            raise BadData(
                'unexpected file names, can not match fwd vs.rev files',
                Err.E16
            )
    else:
        raise BadData('inconsistent number of files / endedness?', Err.E17)

    return [row['Q2'] for _, row in stats.items()]


def basic_checks(data, file_dir):
    """ Checks that apply to any individual summary """
    GOOD = 0.5  # fraction of good alignment

    errors = []
    if data['good_alignments_count'] / data['alignment_winners_count'] < GOOD:
        errors.append((f'{file_dir}: too few good alignments', Err.E1))

    if 'model' not in data:
        errors.append((
            f'{file_dir}: Insufficient summary information: no model',
            Err.E2
        ))
        raise InsufficientData(errors)

    return errors


def check_single(single):
    """
    Consider summaries of single direction

    Returns dictionary with target info.  Raises BadData if it encouters an
    error such that continuing makes no sense.
    """
    MIN_COUNT = 100
    MAX_DISTANCE = 30

    info = {'layout': 'single'}
    errors = []

    info['model_name'] = single['model'].name

    if 'direction' in single:
        # 1. Did reads get swapped?
        if single['direction'] == 'fwd':
            # normal
            info['model_rc'] = False
        elif single['direction'] == 'rev':
            # Reads are RC w.r.t. HMM model
            info['model_rc'] = True
        else:
            raise ValueError('Inconsistent directionality!', Err.E3)

    # 2. basic checks
    if 'fwd_primer' in single:
        info['fwd_primer'] = single['fwd_primer'].name
        if single['fwd_count'] < MIN_COUNT:
            errors.append(('too few good reads with fwd primer', Err.E5))
    else:
        errors.append(('no fwd primer detected', Err.E6))

    if 'rev_primer' in single:
        info['rev_primer'] = single['rev_primer'].name
        if single['rev_count'] < MIN_COUNT:
            errors.append(('too few good reads with rev primer', Err.E7))
    else:
        errors.append(('no rev primer detected', Err.E8))

    if 'fwd_avg_score' in single:
        info['fwd_clean'] = single['fwd_clean']
        info['fwd_distance'] = format_distance(single['fwd_avg_score'])
        if abs(single['fwd_avg_score']['avg']) > MAX_DISTANCE:
            # single read's 5' end must be near a fwd primer
            errors.append(('too far away from fwd primer', Err.E9))

    if 'rev_avg_score' in single:
        info['rev_clean'] = single['rev_clean']
        info['rev_distance'] = format_distance(single['rev_avg_score'])
        if abs(single['rev_avg_score']['avg']) > MAX_DISTANCE:
            # single read's 3' end must be near a rev primer
            errors.append(('too far away from rev primer', Err.E10))

    info['errors'] = errors
    return info


def check_paired(fwd, rev):
    """
    Consider summaries of both direction

    Returns dictionary with target info.  Raises BadData if it encouters an
    error such that continuing makes no sense.
    """
    MIN_COUNT = 100
    MAX_DISTANCE = 30

    info = {'layout': 'paired'}
    errors = []

    # 0. Check if there is agreement on the model
    if fwd['model'] == rev['model']:
        info['model_name'] = fwd['model'].name

        # 1. Did reads get swapped?
        # This check only makes sense if model is not in question
        if fwd['direction'] == 'fwd' and rev['direction'] == 'rev':
            # normal
            swapped = False
        elif fwd['direction'] == 'rev' and rev['direction'] == 'fwd':
            # fwd/rev got swapped -- unswap them right here (!)
            swapped = True
            tmp = fwd
            fwd = rev
            rev = tmp
            del tmp
        else:
            raise BadData('Inconsistent directionality!', Err.E3)
    else:
        errors.append(
            (f'fwd-rev-model-mismatch: {fwd["model"]}/{rev["model"]}', Err.E4)
        )
        swapped = False
    info['swapped'] = swapped

    # These are the forward and reverse primers for the pair:
    fwd_primer = fwd.get('fwd_primer')
    rev_primer = rev.get('rev_primer')

    # 2. basic checks
    if fwd_primer:
        info['fwd_primer'] = fwd_primer.name
        if fwd['fwd_count'] < MIN_COUNT:
            errors.append(('too few good reads with fwd primer', Err.E5))
    else:
        errors.append(('no fwd primer detected', Err.E6))

    if rev_primer:
        info['rev_primer'] = rev_primer.name
        if rev['rev_count'] < MIN_COUNT:
            errors.append(('too few good reads with rev primer', Err.E7))
    else:
        errors.append(('no rev primer detected', Err.E8))

    if 'fwd_avg_score' in fwd:
        info['fwd_clean'] = fwd['fwd_clean']
        info['fwd_distance'] = format_distance(fwd['fwd_avg_score'])
        if abs(fwd['fwd_avg_score']['avg']) > MAX_DISTANCE:
            # fwd read's 5' end must be near a fwd primer
            errors.append(('too far away from fwd primer', Err.E9))

    if 'rev_avg_score' in rev:
        info['rev_clean'] = rev['rev_clean']
        info['rev_distance'] = format_distance(rev['rev_avg_score'])
        if abs(rev['rev_avg_score']['avg']) > MAX_DISTANCE:
            # rev read's 5' end must be near a rev primer
            errors.append(('too far away from rev primer', Err.E10))

    # 4. Calculate clean read positions
    match fwd.get('fwd_clean'):
        case True | None: fwd_start = fwd['hmmfrom_avg']
        case False: fwd_start = fwd_primer.end + 1
        case _: raise ValueError()
    fwd_end = fwd['hmmto_avg']
    if rev_primer:
        if rev_primer.start <= fwd_end:
            info['fwd_rev_readthrough'] = fwd_end - rev_primer.start + 1
            info['fwd_rev_clean'] = False
            fwd_end = rev_primer.start - 1
        else:
            info['fwd_rev_clean'] = True

    match rev.get('rev_clean'):
        case True | None: rev_end = rev['hmmto_avg']
        case False: rev_end = rev_primer.start - 1
        case _: raise ValueError()
    rev_start = rev['hmmfrom_avg']
    if fwd_primer:
        if rev_start <= fwd_primer.end:
            info['rev_fwd_readthrough'] = fwd_primer.end - rev_start + 1
            info['rev_fwd_clean'] = False
            rev_start = fwd_primer.end + 1
        else:
            info['rev_fwd_clean'] = True

    info['fwd_hmmfrom'] = fwd_start
    info['fwd_hmmto'] = fwd_end
    info['rev_hmmfrom'] = rev_start
    info['rev_hmmto'] = rev_end

    # 5. Overlap (or gap?)
    overlap = {
        'start': max(fwd_start, rev_start),
        'end': min(fwd_end, rev_end),
    }
    info['overlap'] = overlap
    overlap_len = overlap['end'] - overlap['start'] + 1
    if overlap_len >= 0:
        info['overlap_pct'] = f'{overlap_len / (rev_end - fwd_start):.0%}'
    else:
        info['gap'] = overlap['start'] - overlap['end']

    info['errors'] = errors
    return info


def format_distance(data):
    avg = data['avg']
    upper = data['plus']
    lower = data['minus']
    sign = '+' if data['avg'] >= 0 else ''
    return f'{sign}{avg:.0f} (+{upper:.0f}/-{lower:.0f})'


if __name__ == '__main__':
    cli()
