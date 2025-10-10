"""
Dispatch dada2 jobs
"""
import argparse
from collections import Counter
from datetime import datetime
import filecmp
from pathlib import Path
import shutil
from tempfile import NamedTemporaryFile


SINGLE_MODE_THRESHOLT = 0.8
""" If there is only a single assigned target (besides UNKNOWN) and
the proportion of samples with this definite target is larger than this
thresholt, then all samples will be automatically assigned to this (mode)
target. """

DEFAULT_OUT_TARGETS = 'amplicon_target_assignments.tsv'
DEFAULT_OUT_SAMPLES = 'sample_files.tsv'

UNKNOWN = 'UNKNOWN'


def main():
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument(
        'project_directory',
        help='The project directory, usually data/projects/<dataset> where the'
             ' last component of the given path is interpreted as the dataset '
             'ID.  This directory is expected to contain an "amplicon" '
             'subdirectory that in turn contains one subdirectory per sample.',
    )
    argp.add_argument(
        'target_tabular_input',
        help='Input file, as made by amplicon.tabular_targets module',
    )
    argp.add_argument(
        '--out-targets',
        metavar='<path>',
        default=DEFAULT_OUT_TARGETS,
        help=f'Path to output sample listing file.  This will a tab-separated '
             f'table listing samples and their amplicon targets.  The file '
             f'should normally be reviewed and amended for correctness and '
             f'completeness.  Defaults to "{DEFAULT_OUT_TARGETS}"',
    )
    argp.add_argument(
        '--out-samples',
        metavar='<path>',
        default='./sample_files.tsv',
        help=f'Path to output sample listing file.  This will be a '
             f'tab-separated table listing samples with their paths to raw '
             f'reads.  Defaults to "{DEFAULT_OUT_SAMPLES}"',
    )
    argp.add_argument(
        '--force',
        action='store_true',
        help='Skip check that raw fastq files exist',
    )
    args = argp.parse_args()

    data = read_input(args.target_tabular_input)
    if not data:
        argp.error('no data in input file?')
    write_assignments(data, args.out_targets)
    write_sample_info(
        data, args.project_directory, args.out_samples, force=args.force,
    )


def read_input(input_file):
    """ Read input from tabular target info input file """
    data = []
    with open(input_file) as ifile:
        header = None
        for line in ifile:
            line = line.rstrip('\n')
            if not line.strip() or line.startswith('#'):
                continue

            if header is None:
                header = line.split('\t')
                continue

            data.append(dict(zip(header, line.split('\t'))))
    return data


def write_assignments(data, output_file):
    """ Partition the samples and write the target assignment output file """
    output_file = Path(output_file)
    stats = Counter()
    for row in data:
        if row['errors']:
            row['target'] = UNKNOWN
        else:
            row['target'] = row['model'] + '.' + row['primers']
        row['override'] = ''
        stats[row['target']] += 1

    stats = stats.most_common()
    print(f'Target stats: {stats}')

    if len(stats) == 2:
        (top_target, top_count), (minor_target, minor_count) = stats
        if minor_target == UNKNOWN:
            if (top_count / len(data)) >= SINGLE_MODE_THRESHOLT:
                print(
                    f'Single-target-mode detected, overriding {minor_count} '
                    f'stray (unknown) target assignments to {top_target}'
                )
            for row in data:
                if row['target'] == minor_target:
                    row['override'] = top_target

    data = sorted(
        data,
        # sort UNKNOWN first
        key=lambda x: '' if x['target'] == UNKNOWN else x['target']
    )

    print(f'Writing {output_file} ... ', end='', flush=True)
    with NamedTemporaryFile(mode='w+') as tmpfile:
        header = ('sample_id', 'target', 'override')
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


def write_sample_info(data, project_dir, output_file, force=False):
    """
    Write the sample info outout file

    This should raise exceptions if any of the paths of fastq files do not
    exist.
    """
    project_dir = Path(project_dir)
    # paths as written to output will be relative to the data root
    glamr_root = project_dir.parent.parent.parent
    outrows = []
    for inrow in data:
        swapped = bool(inrow['swapped'])
        fwd_p_cont, _, rev_p_cont = inrow['primer_content'].partition('/')

        sample_dir = project_dir / 'amplicons' / inrow['sample_id']
        sample_dir = sample_dir.resolve().relative_to(glamr_root.resolve())
        fwd_fq = sample_dir / 'reads' / 'raw_fwd_reads.fastq.gz'
        rev_fq = sample_dir / 'reads' / 'raw_rev_reads.fastq.gz'
        if not (glamr_root / fwd_fq).is_file():
            if not force:
                raise FileNotFoundError(f'no such fastq file: {fwd_fq}')
        if not (glamr_root / rev_fq).is_file():
            if not force:
                raise FileNotFoundError(f'no such fastq file: {rev_fq}')
        outrows.append({
            'sample_id': inrow['sample_id'],
            'sample_dir': sample_dir,
            'fwd_primer_content': fwd_p_cont,
            'rev_primer_content': rev_p_cont,
            'fwd_fastq': rev_fq if swapped else fwd_fq,
            'rev_fastq': fwd_fq if swapped else rev_fq,
        })

    print(f'Writing {output_file} ... ', end='', flush=True)
    with open(output_file, 'w') as ofile:
        cols = list(outrows[0].keys())
        header = '\t'.join(cols) + '\n'
        ofile.write(header)
        for row in outrows:
            ofile.write('\t'.join(str(row[i]) for i in cols))
            ofile.write('\n')
    print('[Done]')


if __name__ == '__main__':
    main()
