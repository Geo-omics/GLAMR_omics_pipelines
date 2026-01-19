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
from .hmm import HMM


SINGLE_MODE_THRESHOLT = 0.8
""" If there is only a single assigned target (besides UNKNOWN) and
the proportion of samples with this definite target is larger than this
thresholt, then all samples will be automatically assigned to this (mode)
target. """

DEFAULT_OUT_TARGETS = 'amplicon_target_assignments.tsv'
DEFAULT_OUT_SAMPLES = 'sample_files.tsv'

UNKNOWN = 'UNKNOWN'
SKIP = 'SKIP'
AUTO_OVERRIDE_PREFIX = '__AUTO__'


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
        '--out-samples',
        metavar='<path>',
        default='./sample_files.tsv',
        help=f'Path to output sample listing file.  This will be a '
             f'tab-separated table listing samples with their paths to raw '
             f'reads.  Defaults to "{DEFAULT_OUT_SAMPLES}"',
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
                    out_samples=args.out_samples,
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
    fastq_files,  # passed by Snakefile, but not used, as we'll get them again
    target_tab,
    project_dir,
    out_assignments=None,
    out_samples=None,
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

    write_assignments(data, out_assignments)
    write_sample_info(get_sample_info(data, project_dir), out_samples)


def load_target_table(input_file):
    """ Read input from tabular target info input file """
    all_data = []
    with open(input_file) as ifile:
        header = None
        for line in ifile:
            line = line.rstrip('\n')
            if not line.strip() or line.startswith('#'):
                continue

            if header is None:
                header = line.split('\t')
                continue

            data = dict(zip(header, line.split('\t')))

            for i in ['swapped', 'fwd_clean', 'rev_clean']:
                if i in data:
                    match data[i]:
                        case 'True': data[i] = True
                        case 'False': data[i] = False
                        case _: raise ValueError('invalid value')

            all_data.append(data)
    return all_data


def write_assignments(data, output_file):
    """ Partition the samples and write the target assignment output file """
    output_file = Path(output_file)
    stats = Counter()
    for row in data:
        if 'errors' in row:
            row['target'] = UNKNOWN
        else:
            row['target'] = '.'.join([
                row['model_name'], row['fwd_primer'], row['rev_primer']
            ])
        row.setdefault('override', '')
        stats[row['target']] += 1

    print('Dispatching...')
    for target, count in stats.most_common():
        print(f'  {count:>4} samples to target {target}')

    if len(stats) == 2:
        (top_target, top_count), (minor_target, minor_count) = stats
        if minor_target == UNKNOWN:
            if (top_count / len(data)) >= SINGLE_MODE_THRESHOLT:
                print(
                    f'Single-target-mode detected, overriding {minor_count} '
                    f'stray (unknown) target assignments to {top_target}'
                )
            for row in data:
                if row['target'] == minor_target and not row['override']:
                    row['override'] = AUTO_OVERRIDE_PREFIX + top_target

    data = sorted(
        data,
        # sort UNKNOWN first
        key=lambda x: '' if x['target'] == UNKNOWN else x['target']
    )

    print(f'Writing {output_file} ... ', end='', flush=True)
    with NamedTemporaryFile(mode='w+') as tmpfile:
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


def get_sample_info(data, project_dir):
    """
    Compile the sample info from the data

        project_dir: A resolved() pathlib.Path.
    """
    # For testing, set env variable to pass the raw reads to dada2
    USE_CLEAN = not bool(environ.get('NO_PRIMER_REMOVAL'))

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
            if row['swapped']:
                # "Fix" the swap
                fwd_infix = 'rev'
                rev_infix = 'fwd'
            else:
                # keep as-is
                fwd_infix = 'fwd'
                rev_infix = 'rev'

            if USE_CLEAN and (row.get('fwd_clean') is False or row.get('fwd_rev_clean') is False):  # noqa:E501
                fwd_fq = f'clean.{fwd_infix}_reads.fastq.gz'
            else:
                fwd_fq = f'raw_{fwd_infix}_reads.fastq.gz'

            if USE_CLEAN and (row.get('rev_clean') is False or row.get('rev_fwd_clean') is False):  # noqa:E501
                rev_fq = f'clean.{rev_infix}_reads.fastq.gz'
            else:
                rev_fq = f'raw_{rev_infix}_reads.fastq.gz'

            info['fwd_fastq'] = samp_dir / 'reads' / fwd_fq
            info['rev_fastq'] = samp_dir / 'reads' / rev_fq

        elif row['layout'] == 'single':
            if USE_CLEAN and (row.get('fwd_clean') is False or row.get('rev_clean') is False):  # noqa:E501
                single_fq = 'clean.single_reads.fastq.gz'
            else:
                single_fq = 'raw_single_reads.fastq.gz'

            info['single_fastq'] = samp_dir / 'reads' / single_fq
        else:
            raise ValueError('invalid layout')

        all_info.append(info)
    return all_info


def write_sample_info(rows, output_file):
    """
    Write the sample info output file

    This should raise exceptions if any of the paths of fastq files do not
    exist.
    """
    print(f'Writing {output_file} ... ', end='', flush=True)
    with open(output_file, 'w') as ofile:
        cols = list(rows[0].keys())
        header = '\t'.join(cols) + '\n'
        ofile.write(header)
        for row in rows:
            ofile.write('\t'.join(str(row[i]) for i in cols))
            ofile.write('\n')
    print('[Done]')


def get_fastq_files(path_template, dataset=None):
    """
    Return list of all fastq files -- helper for some Snakemake input function

    The path_template must accept "dataset" formatting.
    """
    targets_tab = Path(path_template.format(dataset=dataset))
    project_dir = targets_tab.parent.resolve()
    data = load_target_table(targets_tab)
    files = []
    for row in get_sample_info(data, project_dir):
        for key in ['fwd_fastq', 'rev_fastq', 'single_fastq']:
            if key in row:
                files.append(row[key])
    return files


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
            sample_id, target, *override = line.split(maxsplit=3)

            if lnum == 1 and sample_id == 'sample' and target == 'target':
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


def target2dada2_dir(target):
    """
    Get name of dada2 output directory from target identifier string

    This is used in the Snakefile as well as the GLAMR DB loading code.
    """
    hmm, fwdprim, revprim = HMM.parse_target(target)
    return f'dada2.{hmm.format_target_spec(fwdprim, revprim)}'


def spec2targets(spec):
    try:
        print(*HMM.spec2targets(spec), sep='\n')
    except Exception as e:
        raise UsageError(f'Error parsing the given spec {spec=}: {e}') from e


if __name__ == '__main__':
    cli()
