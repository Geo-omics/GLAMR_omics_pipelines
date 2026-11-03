#!/usr/bin/env python3
"""
Guess the amplicon target based on summaries of hmmscan

Output is written to stdout.
"""
import argparse
from collections import Counter
import gzip
from itertools import islice
import json
from pathlib import Path
import re

from ..utils import load_stats, logme, UsageError
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


@logme()
def main(summaries, stats_file, outfile=None, fastqs=None):
    """
    fastqs:
        Optional raw fastq files (fwd, rev), used for targets detected by
        their primer sequences: ITS (no profile HMM) and 12S (MiFish).
    """
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

    errors = []

    if len(data) == 1:
        out_data = {
            'layout': 'single',
            'model_name': data[0]['model'].name
        }
        try:
            errors = basic_checks(*data, 'single')
        except (BadData, InsufficientData) as e:
            errors.append(e.args)
        else:
            try:
                out_data.update(check_single(*data))
            except BadData as e:
                errors.append(e.args)
    elif len(data) == 2:
        out_data = {
            'layout': 'paired',
        }
        check_as_pair = True
        try:
            errors = basic_checks(data[0], 'fwd')
        except (BadData, InsufficientData) as e:
            errors.append(e.args)
            check_as_pair = False
        try:
            errors += basic_checks(data[1], 'rev')
        except (BadData, InsufficientData) as e:
            errors.append(e.args)
            check_as_pair = False
        if check_as_pair:
            try:
                out_data.update(check_paired(*data))
            except BadData as e:
                errors.append(e.args)

        hmm_failed = bool(errors or out_data.get('errors'))
        primer_models = []
        if hmm_failed or looks_like_its(*data):
            # The rRNA HMMs only see the genes flanking ITS, check for ITS
            # primers (both primers at the 5' ends of enough reads)
            primer_models.append(ITS_MODEL)
        # The MiFish primers also amplify bacterial 16S, which can outnumber the
        # fish 12S reads so much that the HMMs see no 12S at all, but all reads
        # start with the 12S primers
        primer_models.append(MITO_12S_MODEL)
        for model_name in primer_models if fastqs else []:
            primer_info = check_primers(model_name, *fastqs, data[0]['length'], data[1]['length'])
            if primer_info is not None:
                out_data = {'layout': 'paired', **primer_info}
                errors = []
                break
    else:
        raise RuntimeError('bug')

    errors += out_data.get('errors', [])
    if errors:
        formatted_errs = []
        for item in errors:
            try:
                msg, err = item
            except (TypeError, ValueError):
                # preserve some flexibility here
                formatted_errs.append(str(item))
            else:
                formatted_errs.append(f'E{err}: {msg}')
        out_data['errors'] = formatted_errs
    elif 'errors' in out_data:
        # don't write out an empty list
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

        for direc, primers in [('fwd', model.fwd_primers), ('rev', model.rev_primers)]:
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

    if 'good_alignments_count' not in data:
        # early return in hmm_summarize.main
        raise InsufficientData(
            f'{file_dir}: Insufficient summary information: there were no '
            f'alignments at all',
            Err.E9
        )

    if 'model' not in data:
        raise InsufficientData(
            f'{file_dir}: Insufficient summary information: no model',
            Err.E2
        )

    errors = []
    if data['good_alignments_count'] / data['alignment_winners_count'] < GOOD:
        errors.append((f'{file_dir}: too few good alignments', Err.E1))

    return errors


def check_single(single):
    """
    Consider summaries of single direction

    Returns dictionary with target info.  Raises BadData if it encouters an
    error such that continuing makes no sense.
    """
    MIN_COUNT = 100
    MAX_DISTANCE = 30

    info = {}
    errors = []

    if 'direction' in single:
        # 1. Which way was the gene sequenced?
        if single['direction'] == 'fwd':
            # normal
            info['model_rc'] = False
        elif single['direction'] == 'rev':
            # Reads are RC w.r.t. HMM model
            info['model_rc'] = True
        else:
            raise BadData('Inconsistent directionality!', Err.E3)

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

    info = {}
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


ITS_MODEL = 'ITS_euk'
MITO_12S_MODEL = '12S_rRNA_mito'
IUPAC = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': '[AG]', 'Y': '[CT]',
    'S': '[CG]', 'W': '[AT]', 'K': '[GT]', 'M': '[AC]', 'B': '[CGT]',
    'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]', 'I': '[ACGT]',
}
COMPLEMENT = str.maketrans('ACGTRYSWKMBDHVNI', 'TGCAYRSWMKVHDBNN')


def looks_like_its(fwd, rev):
    """
    Do the HMM summaries look like ITS reads?

    ITS reads align to the flanking rRNA genes (the ITS anchors), with the fwd
    and rev reads usually on different genes.  Plain 18S or 28S amplicons have
    both directions on the same gene.
    """
    anchors = get_models()[ITS_MODEL].anchors
    models = [getattr(i.get('model'), 'name', None) for i in (fwd, rev)]
    if not all(i in anchors for i in models):
        return False
    return not (models[0] == models[1] and models[0] in ('18S_rRNA_euk', '28S_rRNA_euk'))


def primer_regex(sequence, max_offset):
    """ regex for a primer at a read's 5' end, allowing for a few leading bases """
    return re.compile(
        f'^.{{0,{max_offset}}}' + ''.join(IUPAC[i] for i in sequence.upper())
    )


def read_seqs(fastq, num):
    with gzip.open(fastq, 'rt') as ifile:
        return [line.rstrip('\n') for line in islice(ifile, 1, 4 * num, 4)]


def check_primers(model_name, fwd_fastq, rev_fastq, fwd_len, rev_len):
    """
    Detect a target by its primer sequences at the 5' ends of reads

    Used for ITS, which has no profile HMM, and for 12S, where off-target reads
    can outnumber the target in the HMM summaries.  Returns target info like
    check_paired() does or None if no primer pair of the model is found in
    enough reads.  Coordinates are on the model's frame (for ITS the nominal
    ITS frame).
    """
    NUM_READS = 1000
    MIN_COUNT = 100
    MAX_OFFSET = 10  # leading bases before primer, e.g. heterogeneity spacers

    model = get_models()[model_name]
    reads = {'fwd': read_seqs(fwd_fastq, NUM_READS), 'rev': read_seqs(rev_fastq, NUM_READS)}

    def count(primer, direc):
        pat = primer_regex(primer.sequence, MAX_OFFSET)
        return sum(1 for i in reads[direc] if pat.match(i))

    # best primers and whether the reads are swapped
    fwd_counts = {(p, d): count(p, d) for p in model.fwd_primers for d in reads}
    (fwd_primer, fwd_direc), fwd_count = max(fwd_counts.items(), key=lambda x: x[1])
    swapped = fwd_direc == 'rev'
    rev_direc = 'fwd' if swapped else 'rev'
    rev_counts = {p: count(p, rev_direc) for p in model.rev_primers}
    rev_primer, rev_count = max(rev_counts.items(), key=lambda x: x[1])

    if fwd_count < MIN_COUNT or rev_count < MIN_COUNT:
        print(f'[{model.name}] no primer pair: {fwd_primer.name}: {fwd_count}, '
              f'{rev_primer.name}: {rev_count} of {NUM_READS} reads')
        return None
    print(f'[{model.name}] {fwd_primer.name}: {fwd_count}, {rev_primer.name}: {rev_count} '
          f'of {NUM_READS} reads{" (swapped)" if swapped else ""}')

    # read-through into the other primer (short amplicons)
    def readthrough(primer, direc):
        rc = primer.sequence.upper().translate(COMPLEMENT)[::-1]
        pat = re.compile(''.join(IUPAC[i] for i in rc))
        return sum(1 for i in reads[direc] if pat.search(i))
    fwd_rev_clean = readthrough(rev_primer, fwd_direc) == 0
    rev_fwd_clean = readthrough(fwd_primer, rev_direc) == 0

    if swapped:
        fwd_len, rev_len = rev_len, fwd_len
    fwd_start = fwd_primer.end + 1
    fwd_end = min(fwd_start + fwd_len - len(fwd_primer.sequence) - 1, rev_primer.start - 1)
    rev_end = rev_primer.start - 1
    rev_start = max(rev_end - (rev_len - len(rev_primer.sequence)) + 1, fwd_primer.end + 1)

    info = {
        'model_name': model.name,
        'swapped': swapped,
        'fwd_primer': fwd_primer.name,
        'rev_primer': rev_primer.name,
        'fwd_clean': False,
        'fwd_distance': format_distance({'avg': 0, 'plus': 0, 'minus': 0}),
        'rev_clean': False,
        'rev_distance': format_distance({'avg': 0, 'plus': 0, 'minus': 0}),
        'fwd_rev_clean': fwd_rev_clean,
        'rev_fwd_clean': rev_fwd_clean,
        'fwd_hmmfrom': fwd_start,
        'fwd_hmmto': fwd_end,
        'rev_hmmfrom': rev_start,
        'rev_hmmto': rev_end,
        'overlap': {'start': max(fwd_start, rev_start), 'end': min(fwd_end, rev_end)},
        'errors': [],
    }
    overlap_len = info['overlap']['end'] - info['overlap']['start'] + 1
    if overlap_len >= 0:
        info['overlap_pct'] = f'{overlap_len / (rev_end - fwd_start):.0%}'
    else:
        info['gap'] = info['overlap']['start'] - info['overlap']['end']
    return info


def format_distance(data):
    avg = data['avg']
    upper = data['plus']
    lower = data['minus']
    sign = '+' if data['avg'] >= 0 else ''
    return f'{sign}{avg:.0f} (+{upper:.0f}/-{lower:.0f})'


if __name__ == '__main__':
    cli()
