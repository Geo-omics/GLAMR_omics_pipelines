"""
Summarize results from nhmmscan
"""
import argparse
from collections import Counter
from itertools import groupby
import json
from pathlib import Path
from statistics import median, quantiles

from matplotlib.backends.backend_pdf import PdfPages
import pandas

from ..utils import UsageError
from . import dispatch
from .alignment import HMMRAlignment
from .hmm import HMM


def cli():
    """
    Implement the cmd line tool

    The summary is written to stdout.
    """
    argp = argparse.ArgumentParser(description=__doc__)
    subp = argp.add_subparsers(dest='subcmd')
    sump = subp.add_parser(
        'summary',
        help='Summarize alignments for one sample',
    )
    sump.add_argument(
        'inputfile',
        help='Tabular alignment file as made with nhmmscan\'s --tblout option',
    )
    sump.add_argument(
        '--output',
        help='Name of output file. By default print to stdout.',
    )
    sump.add_argument('--debug', action='store_true', help='debug mode')
    plotp = subp.add_parser(
        'plot',
        help='Plot sequencing depth for a sample',
    )
    plotp.add_argument(
        'inputfile',
        help='Tabular alignment file as made with nhmmscan\'s --tblout option',
    )
    multip = subp.add_parser(
        'multiplot',
        help='Plot sequencing depth for whole dataset',
    )
    multip.add_argument(
        'inputfile',
        help='The target assignment file for the dataset.',
    )
    args = argp.parse_args()
    try:
        match args.subcmd:
            case 'summary': main(args.inputfile, args.output, args.debug)
            case 'plot': plot(args.inputfile)
            case 'multiplot': multiplot(args.inputfile)
            case _: raise ValueError('invalid subcommand')
    except UsageError as e:
        argp.error(e)


def main(infile, outfile=None, debug=False):
    summary = {}

    # 1. load data
    alignments = HMMRAlignment.load(infile)
    summary['total_alignment_count'] = len(alignments)
    if not alignments:
        # nothing else to be done
        finish(summary, outfile)
        return

    # 2. get best alignment per read
    alignments = pick_winning_model(alignments)
    summary['alignment_winners_count'] = len(alignments)
    summary['winners_per_model'] = {}
    for name, count in Counter(i.model.name for i in alignments).most_common():
        summary['winners_per_model'][name] = count

    # 3. filter out "bad" reads
    alignments = [i for i in alignments if i.pass_test]
    summary['good_alignments_count'] = len(alignments)
    if not alignments:
        finish(summary, outfile)
        return

    if debug:
        lines = []
        for i in alignments:
            direc = '->' if i.strand == '+' else '<-'
            out = (i.fwd_match, f'{i.hmmfrom}{direc}{i.hmmto}', i.rev_match)
            lines.append(' '.join(str(i) for i in out))
        for line, count in reversed(Counter(lines).most_common()):
            print(f'{count:>5}', line)
        finish(summary)
        return

    # 4. find primers and mode
    mode = get_mode(alignments)
    top_alignments = mode.pop('top_alignments')
    mode['top_alignment_count'] = len(top_alignments)
    summary.update(mode)
    finish(summary, outfile)


def finish(summary, output=None):
    out_txt = json.dumps(summary, indent=4)
    if output is None:
        print(out_txt)
    else:
        with open(output, 'w') as ofile:
            ofile.write(out_txt)
            ofile.write('\n')
        print(f'Output written to {output}')


def pick_winning_model(alignments):
    """
    Pick best alignment per read

    Usually several model admit some hit, so this picks the best model for each
    read.
    """
    winners = []
    tie_count = 0
    for _, grp in groupby(alignments, key=lambda x: x.qname):
        top = None
        highscore = None
        tiescore = None
        for row in grp:
            if highscore is None or row.score > highscore:
                top = row
                highscore = row.score
            elif row.score == highscore:
                tiescore = row.score
        if tiescore is not None and tiescore == highscore:
            # Top scores are tied, discard this, keeping the earlier alignment.
            # Expecting this to be a rare occurrence, but just in case it's not
            # we're keeping count of ties and report it.
            tie_count += 1
            continue
        winners.append(top)
    if tie_count:
        print(f'# discarding due to tied top hit score: {tie_count}')
    return winners


def get_mode(alignments):
    """ get the most common alignments """
    # Separately get most common fwd and reverse primers, and pairs
    fwd_data = Counter()
    rev_data = Counter()
    pair_data = Counter()
    for i in alignments:
        # below triggers the primer matching
        # (see fwd_match() / rev_match() properties in alignment module)
        fwd_primer = None if i.fwd_match is None else i.fwd_match.primer
        rev_primer = None if i.rev_match is None else i.rev_match.primer
        fwd_data[fwd_primer] += 1
        rev_data[rev_primer] += 1
        pair_data[(fwd_primer, rev_primer, i.direction)] += 1

    top_fwd_primer, fwd_count = fwd_data.most_common()[0]
    top_rev_primer, rev_count = rev_data.most_common()[0]
    top_pair, _ = pair_data.most_common()[0]

    # Check if there is agreement
    if (top_fwd_primer, top_rev_primer) != top_pair[:2]:
        # TODO: what to do?
        print(top_fwd_primer)
        print(top_rev_primer)
        print(top_pair)
        raise NotImplementedError()

    # return data
    items = dict(
        top_alignments=[],
    )

    if top_fwd_primer is not None:
        items['fwd_primer'] = top_fwd_primer.name
        items['fwd_count'] = fwd_count

    if top_rev_primer is not None:
        items['rev_primer'] = top_rev_primer.name
        items['rev_count'] = fwd_count

    clean_fwd_scores = []
    dirty_fwd_scores = []
    clean_rev_scores = []
    dirty_rev_scores = []
    hmmfroms = []
    hmmtos = []
    for i in alignments:
        fwd_primer = None if i.fwd_match is None else i.fwd_match.primer
        rev_primer = None if i.rev_match is None else i.rev_match.primer

        if (fwd_primer, rev_primer, i.direction) == top_pair:
            items['top_alignments'].append(i)
            hmmfroms.append(i.hmmfrom)
            hmmtos.append(i.hmmto)

        if top_fwd_primer is not None and fwd_primer == top_fwd_primer:
            if i.fwd_match.clean:
                clean_fwd_scores.append(i.fwd_match.score)
            else:
                dirty_fwd_scores.append(i.fwd_match.score)

        if top_rev_primer is not None and rev_primer == top_rev_primer:
            if i.rev_match.clean:
                clean_rev_scores.append(i.rev_match.score)
            else:
                dirty_rev_scores.append(i.rev_match.score)

    if len(models := set(i.model for i in items['top_alignments'])) == 1:
        items['model'] = models.pop().name
    else:
        # the pairing logic above should preclude this possibility
        raise RuntimeError('multiple models among top alignments')

    if len(dirs := set(i.direction for i in items['top_alignments'])) == 1:
        items['direction'] = dirs.pop()
    else:
        # the pairing logic above should preclude this possibility
        raise RuntimeError('multiple directions among top alignments')

    if clean_fwd_scores or dirty_fwd_scores:
        items['fwd_clean'] = (len(clean_fwd_scores) > len(dirty_fwd_scores))
        if items['fwd_clean']:
            fwd_scores = clean_fwd_scores
        else:
            fwd_scores = dirty_fwd_scores

        quarts = quantiles(fwd_scores, n=4)
        items['fwd_avg_score'] = {
            'avg': quarts[1],  # median
            'plus': quarts[2] - quarts[1],
            'minus': quarts[1] - quarts[0],
            'count': len(
                [i for i in fwd_scores if quarts[0] <= i <= quarts[2]]
            ),
        }

    if clean_rev_scores or dirty_rev_scores:
        # TODO: review threasholt, when things turn dirty
        items['rev_clean'] = (len(clean_rev_scores) > len(dirty_rev_scores))
        if items['rev_clean']:
            rev_scores = clean_rev_scores
        else:
            rev_scores = dirty_rev_scores

        quarts = quantiles(rev_scores, n=4)
        items['rev_avg_score'] = {
            'avg': quarts[1],  # median
            'plus': quarts[2] - quarts[1],
            'minus': quarts[1] - quarts[0],
            'count': len(
                [i for i in rev_scores if quarts[0] <= i <= quarts[2]]
            ),
        }

    items['hmmfrom_avg'] = median(hmmfroms)
    items['hmmto_avg'] = median(hmmtos)
    return items


def plot(inputfile):
    """ Plot from single HMMR tbl output file """
    alignments = HMMRAlignment.load(inputfile)
    depth = {}
    for a in alignments:
        if a.model not in depth:
            depth[a.model] = {'+': Counter(), '-': Counter()}
        for pos in range(a.hmmfrom, a.hmmto + 1):
            depth[a.model][a.strand][pos] += 1

    for hmm, data in depth.items():
        df_data = {}
        index = pandas.RangeIndex(1, hmm.length + 1)
        for strand, counts in data.items():
            df_data[strand] = pandas.Series(counts, index)
        df = pandas.DataFrame(df_data)
        ax = df.plot()
        ax.set_title(hmm.name)
        ax.set_yscale('log')
        for i in hmm.fwd_primers:
            ax.fill_between(
                [i.start, i.end], 0, 1200,
                color='green', alpha=0.25,
            )
        for i in hmm.rev_primers:
            ax.fill_between(
                [i.start, i.end], 0, 1200,
                color='brown', alpha=0.25,
            )
        ax.figure.savefig(f'{hmm.name}.depth.pdf')
        ax.cla()


def multiplot(inputfile):
    """ Plot from multiple HMMR tbl output files """
    infile = Path(inputfile)
    samp_dir = infile.parent / 'amplicons'
    groupings = None
    if infile.is_file():
        assignments = dispatch.get_assignments(path=infile)
        groupings = {}
        for sample_id, target in assignments.items():
            if target in (dispatch.SKIP, dispatch.UNKNOWN):
                groupings[sample_id] = target
            else:
                hmm, _, _ = HMM.parse_target(target)
                groupings[sample_id] = hmm
    elif samp_dir.is_dir():
        print('[NOTICE] no such assignmetn file entering dataset-auto-mode')
    else:
        raise UsageError(f'no such file {inputfile} and no such directory '
                         f'{samp_dir}')

    data = {}
    for i in samp_dir.glob('samp_*/detect_region/tbl_single.txt'):
        sample_id = i.parts[-3]
        if groupings is None:
            group = 'all'
        else:
            group = groupings[sample_id]

        if group not in data:
            data[group] = []

        grpdat = {}
        for a in HMMRAlignment.load(i):
            if isinstance(group, str):
                # keep data for any and all models
                pass
            else:
                # group is a HMM
                if a.model != group:
                    # skip other model's alignments
                    continue
            if a.model not in grpdat:
                grpdat[a.model] = []
            depth = {'+': Counter(), '-': Counter()}
            for pos in range(a.hmmfrom, a.hmmto + 1):
                depth[a.strand][pos] += 1
            grpdat[a.model].append((sample_id, depth))
        data[group] = grpdat

    pp = PdfPages('.pdf')
    for group, grpdat in data.items():
        for hmm, depths in grpdat.items():
            index = pandas.RangeIndex(1, hmm.length + 1)

            df_data = {}
            for sample_id, depth in depths:
                for strand, counts in depth.items():
                    df_data[f'{sample_id}{strand}'] = \
                        pandas.Series(counts, index)

            df = pandas.DataFrame(df_data)
            ax = df.plot()
            ax.set_title(hmm.name)
            ax.set_yscale('log')
            for i in hmm.fwd_primers:
                ax.fill_between(
                    [i.start, i.end], 0, 1200,
                    color='green', alpha=0.25,
                )
            for i in hmm.rev_primers:
                ax.fill_between(
                    [i.start, i.end], 0, 1200,
                    color='brown', alpha=0.25,
                )
            ax.figure.savefig(f'{hmm.name}.depth.pdf')
            ax.cla()



if __name__ == '__main__':
    cli()
