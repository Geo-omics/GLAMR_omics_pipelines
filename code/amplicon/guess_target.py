#!/usr/bin/env python3
"""
Guess the amplicon target based on summaries of hmmscan

Output is written to stdout.
"""
import argparse
import json

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

    args = argp.parse_args()
    try:
        main(args.summary, args.stats, args.output)
    except UsageError as e:
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
        lengths = get_lengths(stats_file)
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
    out_txt = json.dumps(out_data, indent=4)
    if outfile is None:
        print(out_txt)
    else:
        with open(outfile, 'w') as ofile:
            ofile.write(out_txt)
            ofile.write('\n')
        print(f'Output written to: {outfile}')


class UsageError(Exception):
    pass


class BadData(Exception):
    pass


class InsufficientData(Exception):
    pass


def error(errno, msg):
    print(f'error E{errno.value}:', msg)


def read_summary(file):
    with open(file) as ifile:
        return json.load(ifile)


def get_lengths(stats_file):
    """
    Get length(s) from stats file

    Assume first row is forward data and second row i(if any) reverse.  Get the
    median (truncated towards zero) as nominal read length.
    """
    lengths = []
    with open(stats_file) as ifile:
        colidx = None
        for line in ifile:
            row = line.rstrip('\n').split('\t')
            if colidx is None:
                # this is the header, first row
                try:
                    colidx = row.index('Q2')
                except ValueError as e:
                    raise BadData('failed parsing header', Err.E15) from e

                continue

            try:
                lengths.append(int(float(row[colidx])))
            except (ValueError, TypeError, KeyError) as e:
                raise BadData(
                    f'Failed reading/parsing stats file: '
                    f'{e.__class__.__name__}: {e}',
                    Err.E16
                ) from e

    if len(lengths) == 0:
        raise BadData('stats data row missing', Err.E17)
    if len(lengths) > 2:
        raise BadData('too many stats data rows', Err.E18)

    return lengths


def basic_checks(data, file_dir):
    """ Checks that apply to any individual summary """
    GOOD = 0.5  # fraction of good alignment

    errors = []
    if data['good_alignments_count'] / data['alignment_winners_count'] < GOOD:
        errors.append(f'{file_dir}: too few good alignments', Err.E1)

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
    models = get_models()
    errors = []

    info['model_name'] = models[single['model']].name

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
        info['fwd_primer'] = single['fwd_primer']
        if single['fwd_count'] < MIN_COUNT:
            errors.append(('too few good reads with fwd primer', Err.E5))
    else:
        errors.append(('no fwd primer detected', Err.E6))

    if 'rev_primer' in single:
        info['rev_primer'] = single['rev_primer']
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
    models = get_models()
    errors = []

    # 0. Check if there is agreement on the model
    if fwd['model'] == rev['model']:
        info['model_name'] = models[fwd['model']].name

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

    # 2. basic checks
    if 'fwd_primer' in fwd:
        info['fwd_primer'] = fwd['fwd_primer']
        if fwd['fwd_count'] < MIN_COUNT:
            errors.append(('too few good reads with fwd primer', Err.E5))
    else:
        errors.append(('no fwd primer detected', Err.E6))

    if 'rev_primer' in rev:
        info['rev_primer'] = rev['rev_primer']
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

    # 4. Calculate overlap
    if 'fwd_avg_scores' in fwd and 'rev_avg_scores' in rev:
        if fwd['fwd_clean']:
            start = fwd['fwd_primer'].end + fwd['fwd_avg_score']['avg']
        else:
            start = fwd['fwd_primer'].end
        if rev['rev_clean']:
            end = rev['rev_primer'].start + rev['rev_avg_score']['avg']
        else:
            end = rev['rev_primer'].start
        overlap = start + fwd['length'] + rev['length'] - end
        info['overlap'] = overlap
        if overlap >= 0:
            info['overlap_pct'] = f'{overlap / (end - start):.0%}'

    else:
        overlap = None

    if overlap is not None and overlap >= 0:
        # check running into primer at other ends
        if 'rev_primer' in fwd and fwd['rev_primer'] == rev['rev_primer']:
            if 'rev_clean' in fwd and not fwd['rev_clean']:
                info['fwd_reads_run_into_rev_primer'] = fwd['rev_avg_score']['avg']  # noqa: E501
        if 'fwd_primer' in rev and rev['fwd_primer'] == fwd['fwd_primer']:
            if 'fwd_clean' in rev and not rev['fwd_clean']:
                info['rev_reads_run_into_fwd_primer'] = rev['fwd_avg_score']['avg']  # noqa: E501

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
