"""
Dispatch dada2 jobs for a dataset
"""
import argparse
from collections import Counter
from datetime import datetime
import filecmp
from os import environ
from pathlib import Path
import shutil
from tempfile import NamedTemporaryFile

from ..utils import UsageError
from . import get_models
from .hmm import HMM


SINGLE_MODE_THRESHOLT = 0.8
""" If there is only a single assigned target (besides UNKNOWN) and
the proportion of samples with this definite target is larger than this
thresholt, then all samples will be automatically assigned to this (mode)
target. """

DEFAULT_OUT_TARGETS = 'amplicon_target_assignments.tsv'
DEFAULT_OUT_DADA2_INPUT = 'dada2_input.tsv'

UNKNOWN = 'UNKNOWN'
SKIP = 'SKIP'
AUTO_OVERRIDE_PREFIX = '__AUTO__'
PRIMER_SEP = '<->'


def cli():
    argp = argparse.ArgumentParser(description=__doc__)
    subp = argp.add_subparsers(
        dest='subcommand',
    )
    makep = subp.add_parser(
        'make',
        help='Do the dispatching, so make the three dispatch output files.',
    )
    makep.add_argument(
        'project_directory',
        help='The project directory, usually data/projects/<dataset> where the'
             ' last component of the given path is interpreted as the dataset '
             'ID.  This directory is expected to contain an "amplicon" '
             'subdirectory that in turn contains one subdirectory per sample.',
    )
    makep.add_argument(
        'target_tabular_input',
        help='Input file, as made by amplicon.tabular_targets module',
    )
    makep.add_argument(
        '--out-assignments',
        metavar='<path>',
        default=DEFAULT_OUT_TARGETS,
        help=f'Path to target assignment output file.  This will a tab-'
             f'separated table listing samples and their amplicon targets.  '
             f'The file should normally be reviewed and amended for '
             f'correctness and completeness.  Defaults to '
             f'"{DEFAULT_OUT_TARGETS}"',
    )
    makep.add_argument(
        '--out-dada2-input',
        metavar='<path>',
        default='./dada2_input.tsv',
        help=f'Path to output dada2-input file.  This will be a '
             f'tab-separated table listing samples with their paths to raw '
             f'reads and other bits.  Defaults to "{DEFAULT_OUT_DADA2_INPUT}"',
    )
    spec2targetsp = subp.add_parser(
        'spec2targets',
        help='List possible primer combinations for given target spec',
    )
    spec2targetsp.add_argument('spec')
    argp.add_argument('--traceback', action='store_true',
                      help='always print traceback on error')
    args = argp.parse_args()
    try:
        match args.subcommand:
            case 'make':
                make(
                    None,
                    args.target_tabular_input,
                    args.project_directory,
                    out_assignments=args.out_assignments,
                    out_dada2_input=args.out_dada2_input,
                )
            case 'spec2targets':
                spec2targets(args.spec)
            case _: raise ValueError('invalid subcommand')
    except UsageError as e:
        if args.traceback:
            raise e
        else:
            argp.error(e)


def make(
    target_tab,
    project_dir,
    out_assignments=None,
    out_dada2_input=None,
):
    project_dir = Path(project_dir).resolve()
    data = load_target_table(target_tab)
    if not data:
        raise UsageError(f'No data in input file: {target_tab}')

    # sort by sample ID
    def by_samp_id_key(row):
        value = row['sample']
        try:
            value = int(value.removeprefix('samp_'))
        finally:
            return value
    data = sorted(data, key=by_samp_id_key)

    if (out_assignments := Path(out_assignments)).is_file():
        # honor existing SKIP overrides
        overrides = get_assignments(
            path=out_assignments,
            manual_overrides_only=True,
        )
        for row in data:
            if row['sample'] in overrides:
                row['override'] = overrides[row['sample']]

    groups = dispatch(data)
    write_assignments(data, groups, out_assignments)
    apply_overrides(data)
    write_dada2_input(get_dada2_input(data, project_dir), out_dada2_input)


def load_target_table(input_file):
    """ Read input from tabular target info input file """
    hmms = get_models()
    rows = []
    with open(input_file) as ifile:
        header = None
        for lnum, line in enumerate(ifile, start=1):
            line = line.rstrip('\n')
            if not line.strip() or line.startswith('#'):
                continue

            if header is None:
                header = line.split('\t')
                continue

            row = dict(zip(header, line.split('\t')))

            # Convert row data in certain boolean columns
            for i in ['swapped', 'fwd_clean', 'rev_clean']:
                if i in row:
                    match row[i]:
                        case 'True': row[i] = True
                        case 'False': row[i] = False
                        case _: row[i] = None

            if row.get('errors'):
                row['target'] = UNKNOWN
            else:
                hmm = hmms[row['model_name']]
                row['target'] = hmm.format_target_spec(
                    hmm.fwd_primer[row['fwd_primer']],
                    hmm.rev_primer[row['rev_primer']]
                )
            row.setdefault('override', '')
            row.setdefault('is_auto_override', False)
            rows.append(row)
    return rows


def dispatch(data):
    """ Partition the samples and apply auto overrides """
    groups = Counter()
    for row in data:
        if row.get('errors'):
            groups[(UNKNOWN, None, None)] += 1
        else:
            groups[(row['target'], row['fwd_primer'], row['rev_primer'])] += 1

    print('Dispatching...')
    for (target, fwd_p, rev_p), count in groups.most_common():
        if target is UNKNOWN:
            target_str = UNKNOWN
        else:
            target_str = f'{target}:  {fwd_p}{PRIMER_SEP}{rev_p}'
        print(f'  {count:>4} samples to target {target_str}')

    if len(groups) == 2:
        (top_target, top_count), (minor_target, minor_count) = groups.items()
        if minor_target == (UNKNOWN, None, None):
            if (top_count / len(data)) >= SINGLE_MODE_THRESHOLT:
                print(
                    f'Single-target-mode detected, overriding {minor_count} '
                    f'stray (unknown) target assignments to {top_target}'
                )
            top_target, top_fwd_p, top_rev_p = top_target
            for row in data:
                if row['target'] == UNKNOWN and not row['override']:
                    row['override'] = f'{top_target}:{top_fwd_p}<>{top_rev_p}'
                    row['is_auto_override'] = True
    return groups


def write_assignments(data, groups, output_file):
    """ Write the target assignment file """
    for row in data:
        if row['is_auto_override']:
            row['override'] = AUTO_OVERRIDE_PREFIX + row['override']

    data = sorted(
        data,
        # sort UNKNOWN first
        key=lambda x: '' if x['target'] == UNKNOWN else x['target']
    )

    output_file = Path(output_file)
    print(f'Writing {output_file} ... ', end='', flush=True)
    with NamedTemporaryFile(mode='w+') as tmpfile:
        tmpfile.write(
            '# *** This file assigns each sample to an amplicon target ***\n'
            '# How to manually override: write a target spec, with or without\n'
            '# primers (which, if specified, are scheduled to be removed) into the\n'
            '# third column.\n'
            '# The format is: HMM_NAME.START-END:FWD_PRIMER_NAME<->REV_PRIMER_NAME\n'
            '# The primer spec, starting w/colon is optional.  Each primer may be\n'
            '# a single primer name or a comma-separated list of such primer names.\n'
            '# You may add additional comment lines as needed.  After editing this\n'
            '# file, you must delete the dada2_input file and any dada2 output for\n'
            '# dada2 runs that would be affected by your changes and then re-run\n'
            '# the workflow.\n'
        )
        example_lines = []
        for (target, fwd_p, rev_p) in groups:
            if target == UNKNOWN:
                continue
            example_lines.append(f'#    {target}:{fwd_p}{PRIMER_SEP}{rev_p}\n')
        if example_lines:
            tmpfile.write('# Detected targets:\n')
            tmpfile.writelines(example_lines)

        header = ('sample', 'target', 'override')
        tmpfile.write('\t'.join(header))
        tmpfile.write('\n')
        for row in data:
            tmpfile.write('\t'.join((row[i] for i in header)))
            tmpfile.write('\n')
        tmpfile.seek(0)
        if output_file.exists():
            if filecmp.cmp(output_file, tmpfile.name, shallow=False):
                print('already is up-to-date [OK]')
                return
            else:
                # make backup copy
                mtime = (
                    datetime
                    .fromtimestamp(output_file.stat().st_mtime)
                    .isoformat()
                )
                backup = output_file.with_name(output_file.name + '.' + mtime)
                print(f'making backup: {backup} ... ', end='', flush=True)
                output_file.rename(backup)
        with open(output_file, 'w') as ofile:
            tmpfile.seek(0)
            shutil.copyfileobj(tmpfile, ofile)

    print('[Done]')


def apply_overrides(data):
    """ Helper to change bits according to overrides """
    for row in data:
        if not row['override']:
            continue

        target, _, primers = row['override'].partition(':')
        fwd_p, _, rev_p = primers.partition(PRIMER_SEP)
        row['target'] = target.strip()
        row['fwd_clean'] = True if fwd_p.strip() == '' else False
        row['rev_clean'] = True if rev_p.strip() == '' else False
        row['fwd_rev_clean'] = True  # FIXME need more general solution for read-through
        row['rev_fwd_clean'] = True  # FIXME need more general solution for read-through


def get_dada2_input(data, project_dir):
    """
    Compile the dada input from the data

        project_dir: A resolved() pathlib.Path.
    """
    # For testing, set env variable to pass the raw reads to dada2
    USE_CLEAN = not bool(environ.get('NO_PRIMER_REMOVAL'))

    EXTRA_KEYS_PAIRED = [
        'fwd_hmmfrom', 'fwd_hmmto', 'rev_hmmfrom', 'rev_hmmto', 'overlap_start',
        'overlap_end',
    ]
    """ bits to be copied, as needed by code/ampliconTrunc.R """

    # paths as written to output will be relative to the data root
    glamr_root = project_dir.parent.parent.parent
    all_info = []
    for row in data:
        samp_dir = project_dir / 'amplicons' / row['sample']
        samp_dir = samp_dir.resolve().relative_to(glamr_root)

        info = {
            'sample': row['sample'],
            'sample_dir': samp_dir,
        }

        if row['layout'] == 'paired':
            # set default blank values
            info['fwd_fastq'] = None
            info['rev_fastq'] = None
            for key in EXTRA_KEYS_PAIRED:
                info[key] = None

            if None not in (row['swapped'], row['fwd_clean'], row['rev_clean']):
                if row['swapped']:
                    # "Fix" the swap
                    fwd_infix = 'rev'
                    rev_infix = 'fwd'
                    if not row['override']:
                        for key in EXTRA_KEYS_PAIRED:
                            if 'fwd' in key:
                                key2 = key.replace('fwd', 'rev')
                            elif 'rev' in key:
                                key2 = key.replace('rev', 'fwd')
                            else:
                                key2 = key
                            info[key2] = row[key]
                else:
                    # keep as-is
                    fwd_infix = 'fwd'
                    rev_infix = 'rev'
                    if not row['override']:
                        for key in EXTRA_KEYS_PAIRED:
                            info[key] = row[key]

                # remove_primers already un-swaps the clean reads, so only the
                # raw reads get the swap fix
                if USE_CLEAN and (row.get('fwd_clean') is False or row.get('fwd_rev_clean') is False):  # noqa:E501
                    fwd_fq = 'clean.fwd_reads.fastq.gz'
                else:
                    fwd_fq = f'raw_{fwd_infix}_reads.fastq.gz'

                if USE_CLEAN and (row.get('rev_clean') is False or row.get('rev_fwd_clean') is False):  # noqa:E501
                    rev_fq = 'clean.rev_reads.fastq.gz'
                else:
                    rev_fq = f'raw_{rev_infix}_reads.fastq.gz'

                info['fwd_fastq'] = samp_dir / 'reads' / fwd_fq
                info['rev_fastq'] = samp_dir / 'reads' / rev_fq

        elif row['layout'] == 'single':
            if None in (row['fwd_clean'], row['rev_clean']):
                info['single_fastq'] = None
            else:
                if USE_CLEAN and (row.get('fwd_clean') is False or row.get('rev_clean') is False):  # noqa:E501
                    single_fq = 'clean.single_reads.fastq.gz'
                else:
                    single_fq = 'raw_single_reads.fastq.gz'

                info['single_fastq'] = samp_dir / 'reads' / single_fq
        else:
            raise ValueError('invalid layout')

        all_info.append(info)
    return all_info


def write_dada2_input(rows, output_file):
    """ Write the dada2-input file """
    print(f'Writing {output_file} ... ', end='', flush=True)
    with open(output_file, 'w') as ofile:
        cols = list(rows[0].keys())
        header = '\t'.join(cols) + '\n'
        ofile.write(header)
        for row in rows:
            ofile.write('\t'.join(str(row[i] or '') for i in cols))
            ofile.write('\n')
    print('[Done]')


def get_fastq_files(dada2_input_tmpl, assignments_tmpl, **wildcards):
    """
    Return list of all fastq files for a dada2 run

    Helper for some Snakemake input function

    The path templates must accept the wildcard formattings as needed.  Returns
    a list of paths, as found in the fastq column(s).
    """
    assignments = get_assignments(path=Path(assignments_tmpl.format(**wildcards)))
    dada2_input = Path(dada2_input_tmpl.format(**wildcards))
    with open(dada2_input) as ifile:
        head = ifile.readline().strip().split('\t')
        if head[0] != 'sample':
            raise RuntimeError('first column must be "sample"')

        # figure out where the fastq paths are stored
        col_idxs = []
        try:
            col_idxs.append(head.index('fwd_fastq'))
            col_idxs.append(head.index('rev_fastq'))
        except ValueError as e1:
            try:
                col_idxs.append(head.index('single_fastq'))
            except ValueError as e2:
                raise RuntimeError(f'Unexpected file header: {e1} and {e2}') from e1

        file_list = []
        for line in ifile:
            row = line.strip().split('\t')
            if target_spec := wildcards.get('target_spec'):
                sample = row[0]
                if target_spec != assignments[sample].split(':', maxsplit=1)[0]:
                    # target does not match
                    continue
            for i in col_idxs:
                file_list.append(row[i])

    return file_list


def get_assignments(
    dataset_id=None,
    omics_pipeline_root=None,
    path=None,
    manual_overrides_only=False,
):
    """ Parse the assignments file """
    if dataset_id is None and path is None:
        raise ValueError('dataset_id or path is required')

    if dataset_id:
        if omics_pipeline_root is None:
            omics_pipeline_root = Path()
        path = (omics_pipeline_root / 'data' / 'projects' / dataset_id
                / DEFAULT_OUT_TARGETS)
    else:
        # path is given
        pass

    data = {}
    with open(path) as ifile:
        for lnum, line in enumerate(ifile, start=1):
            # The parsing below tolerates a bit, since the file may be manually
            # edited
            line = line.strip()
            if line.startswith('#') or not line:
                continue

            # require two fields, override is an optional third field,
            # ignore anything beyond that (could be used for a comment)
            # Delimiter is any blank space
            sample_id, target, *override = line.split(maxsplit=3)

            if not data and sample_id == 'sample' and target == 'target':
                # header line
                continue

            if manual_overrides_only:
                if override:
                    if not override[0].startswith(AUTO_OVERRIDE_PREFIX):
                        data[sample_id] = ' '.join(override)  # incl. comments
                continue

            override = override[0] if override else ''

            target = override.removeprefix(AUTO_OVERRIDE_PREFIX) or target

            if not target:
                raise RuntimeError(f'empty target at line {lnum} in {path}')

            if sample_id in data:
                raise RuntimeError('duplicate sample at line {lnum} in {path}')

            data[sample_id] = target
    return data


def target2spec(target):
    """
    Get formatted spec (HMM+range) from target (HMM+primers) identifier string

    This is used in the Snakefile as well as the GLAMR DB loading code.
    """
    hmm, fwdprim, revprim = HMM.parse_target(target)
    return hmm.format_target_spec(fwdprim, revprim)


def spec2targets(spec):
    try:
        hmm, fwd_primers, rev_primers = HMM.spec2targets(spec)
    except Exception as e:
        raise UsageError(f'Error parsing the given spec {spec=}: {e}') from e

    if not fwd_primers and not rev_primers:
        print(hmm.name)
    else:
        for i in fwd_primers:
            for j in rev_primers:
                print(f'{hmm.name}.{i.name}.{j.name}')


if __name__ == '__main__':
    cli()
