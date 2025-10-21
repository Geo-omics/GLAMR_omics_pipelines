#!/usr/bin/env python3
"""
Guess the amplicon target based on summaries of hmmscan

Output is written to stdout.
"""
import argparse

from amplicon import get_models, Err


def main():
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument(
        'summary',
        nargs='+',
        type=argparse.FileType(),
        help='One or more paths to summarized nhmmscan results, as made by the'
             ' amplicon_hmm_summarize script.  For paired-end sequencing two '
             'files should be provided, the first being for the forward reads '
             'and the second for the reverse reads.',
    )
    argp.add_argument(
        '--stats',
        type=argparse.FileType(),
        help='fastq statistics file as made with "seqkit stats -a" from raw '
             'fastq files',
    )

    args = argp.parse_args()

    try:
        data = [read_summary(i) for i in args.summary]
    except BadData as e:
        argp.error(f'Error parsing input data: {e}')

    if len(data) == 1:
        raise NotImplementedError('single-ended case not implemented')
    elif len(data) == 2:
        try:
            data[0]['length'], data[1]['length'] = get_lengths(args.stats)
        except BadData as e:
            msg, errno = e.args
            error(errno, f'failed getting read statistics: {msg}')
            argp.exit()

        try:
            check_paired(*data)
        except InsufficientData as e:
            msg, errno = e.args
            error(errno, f'{msg}')
            if fwd_comments := ' // '.join(data[0]['comments']):
                print(f'fwd_summary_extra_info: {fwd_comments}')
            if rev_comments := ' // '.join(data[1]['comments']):
                print(f'rev_summary_extra_info: {rev_comments}')
            argp.exit()
    else:
        args.error('too many arguments, provide at more two summaries')


class BadData(Exception):
    pass


class InsufficientData(Exception):
    pass


def error(errno, msg):
    print(f'error E{errno.value}:', msg)


def read_summary(file):
    data = {}
    comments = []

    for line in file:
        line = line.strip()
        if not line:
            continue
        if line.startswith('#'):
            comments.append(line.lstrip('#').strip())
            continue
        key, _, value = line.partition(':')
        key = key.strip()
        value = value.strip()
        if key in data:
            raise BadData(f'duplicate key: {key}', Err.E12)
        else:
            if not key.isidentifier():
                raise BadData(f'invalid key: {key}', Err.E13)
        if value:
            try:
                value = int(value)
            except ValueError:
                pass
        else:
            value = None
        data[key] = value

    data['comments'] = comments
    return data


def get_lengths(stats_file):
    """
    Get length from stats file

    Assume first row is forward data and second row reverse.  Get the
    median (truncated towards zero) as nominal read length.
    """
    try:
        header, fwdline, revline = (line.rstrip('\n') for line in stats_file)
    except Exception as e:
        raise BadData(
            f'Failed reading/parsing stats file: {e.__class__.__name__}: {e}',
            Err.E14
        ) from e
    try:
        colidx = header.split('\t').index('Q2')
    except ValueError as e:
        raise BadData('failed parsing header', Err.E16) from e

    try:
        fwdlen = int(float(fwdline.split('\t')[colidx]))
        revlen = int(float(revline.split('\t')[colidx]))
    except Exception as e:
        raise BadData(
            f'Failed parsing length data: {e.__class__.__name__}: {e}',
            Err.E17
        ) from e
    return fwdlen, revlen


REQUIRED_ITEMS = (
    'model', 'direction', 'count', 'primer', 'primer_content', 'distance'
)


def check_paired(fwd, rev):
    MIN_COUNT = 100
    MAX_DISTANCE = 5

    missing_fwd = missing_rev = []
    for i in REQUIRED_ITEMS:
        missing_fwd = [i for i in REQUIRED_ITEMS if i not in fwd]
        missing_rev = [i for i in REQUIRED_ITEMS if i not in rev]
    if missing_fwd or missing_rev:
        raise InsufficientData(
            f'Insufficient summary information: Missing these: fwd: '
            f'{missing_fwd} / rev: {missing_rev}',
            Err.E4
        )

    models = get_models()
    errors = []

    # 0. Check if there is agreement on the model
    if fwd['model'] == rev['model']:
        model_name = models[fwd['model']].name

        # 1. Did reads get swapped?
        # This check only makes sense if model is not in question
        if fwd['direction'] == 'fwd' and rev['direction'] == 'rev':
            # normal
            swapped = False
        elif fwd['direction'] == 'rev' and rev['direction'] == 'fwd':
            # fwd/rev got swapped -- unswap them right here
            swapped = True
            tmp = fwd
            fwd = rev
            rev = tmp
            del tmp
        else:
            raise BadData('Inconsistent directionality!', Err.E9)
    else:
        errors.append(('fwd-rev-model-mismatch', Err.E1))
        model_name = f'{fwd["model"]}/{rev["model"]} ???'
        swapped = False

    # 2. basic checks
    if fwd['count'] < MIN_COUNT or rev['count'] < MIN_COUNT:
        errors.append(('too few good alignments', Err.E2))
    if abs(fwd['distance']) > MAX_DISTANCE or abs(rev['distance']) > MAX_DISTANCE:  # noqa:E501
        errors.append(('too far away from primer(s)', Err.E3))

    # 3. get and check primers: consistent with direction?
    for prim in models[fwd['model']].fwd_primers:
        if prim.name == fwd['primer']:
            fwd_primer = dict(name=prim.name, start=prim.start, end=prim.end)
            break
    else:
        raise RuntimeError(f'unknown fwd primer: {fwd["primer"]}')

    for prim in models[rev['model']].rev_primers:
        if prim.name == rev['primer']:
            rev_primer = dict(name=prim.name, start=prim.start, end=prim.end)
            break
    else:
        raise RuntimeError(f'unknown rev primer: {rev["primer"]}')

    # 4. Calculate overlap
    if fwd['model'] == rev['model']:
        if fwd['primer_content']:
            start = fwd_primer['start'] + fwd['distance']
        else:
            start = fwd_primer['end'] + fwd['distance']
        if rev['primer_content']:
            end = rev_primer['end'] - rev['distance']  # TODO: check math
        else:
            end = rev_primer['end'] + rev['distance']  # TODO: check math
        overlap = start + fwd['length'] + rev['length'] - end
        if overlap >= 0:
            overlap_pct = f'{overlap / (end - start):.0%}'
        else:
            overlap_pct = None
    else:
        overlap = overlap_pct = None

    # 5. Compile output
    if fwd['primer_content'] or rev['primer_content']:
        has_fwd_primer = 'YES' if fwd['primer_content'] == 'yes' else 'no'
        has_rev_primer = 'YES' if rev['primer_content'] == 'yes' else 'no'
        primer_content = f'{has_fwd_primer}/{has_rev_primer}'
    else:
        primer_content = None
    if fwd['length'] == rev['length']:
        length = fwd['length']
    else:
        length = f'{fwd["length"]}/{rev["length"]}'

    data = {
        'primers': f'{fwd_primer["name"]}/{rev_primer["name"]}',
        'swapped': 'YES' if swapped else None,
        'model': model_name,
        'primer_content': primer_content,
        'distance': f'{fwd["distance"]}/{rev["distance"]}',
        'length': length,
        'overlap': overlap,
        'overlap_pct': overlap_pct,
    }
    for key, value in data.items():
        if value is None:
            continue
        print(f'{key}: {value}')
    for msg, errno in errors:
        error(errno, msg)


if __name__ == '__main__':
    main()
