"""
Collect targeting information for all samples of a dataset
"""
import argparse
from contextlib import ExitStack
import json
from pathlib import Path
from statistics import mean
import sys

from ..utils import logme, UsageError


TEST_FILE_PATTERN = 'guess*{sample_id}.txt'


def cli():
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument(
        'project_directory',
        help='The project directory, usually data/projects/<dataset> where the'
             ' last component of the given path is interpreted as the dataset '
             'ID.  This directory is expected to contain an "amplicons" '
             'subdirectory that in turn contains one subdirectory (or symlink)'
             ' per sample.',
    )
    argp.add_argument(
        '--outfile',
        help='Path to output file',
    )
    args = argp.parse_args()
    try:
        main(
            None,
            output=args.outfile,
            project_dir=args.project_directory,
        )
    except UsageError as e:
        argp.error(e)


INFILE_TAIL = Path('detect_region') / 'target_info.json'


@logme()
def main(infiles, output=None, dataset=None, project_dir=None):
    project_dir = Path(project_dir)
    if not project_dir.is_dir():
        raise UsageError(f'no such directory: {project_dir}')

    # 1. get input files if needed
    if infiles is None:
        infiles = find_info_files(project_dir)
    else:
        infiles = [Path(i) for i in infiles]

    # extract sample IDs if possible
    samp_ids = []
    tail_len = len(INFILE_TAIL.parts)
    for i in infiles:
        if i.parts[-tail_len:] == INFILE_TAIL.parts:
            samp_ids.append(i.parents[tail_len - 1].name)
        else:
            samp_ids.append(i)

    data = zip(samp_ids, collect_target_data(infiles), strict=True)

    # sort data by target and numeric sample ID
    def sort_key(item):
        samp_id, row = item
        try:
            snum = int(samp_id.removeprefix('samp_'))
        finally:
            snum = samp_id
        return (row['layout'], row['model_name'], row['fwd_primer'],
                row['rev_primer'], snum)
    data = sorted(data, key=sort_key)

    write_table(data, get_columns(data), output=output)
    return
    agg_data = aggregate(data, dataset=dataset)
    write_output(agg_data, output)


def find_info_files(project_dir):
    """ Get all target_info files for a given dataset """
    infiles = []
    for i in (project_dir / 'amplicons').iterdir():
        if not i.is_dir():
            continue

        if i.name.startswith('samp_'):
            if (file := i / INFILE_TAIL).is_file():
                infiles.append(file)
            else:
                print(f'[WARNING] target info missing: {file}')
    return infiles


def collect_target_data(infiles):
    """
    Collect data from target_info files
    """
    data = []
    for i in infiles:
        with open(i) as ifile:
            data.append(flatten(json.load(ifile)))
    return data


def flatten(data):
    """ flatten dict-of-dicts into shallow dict """
    ret = {}
    for k, v in data.items():
        if isinstance(v, dict):
            for kk, vv in v.items():
                ret[f'{k}_{kk}'] = vv
        else:
            ret[k] = v
    return ret


def get_columns(data):
    ranks = {}
    for _, sample_data in data:
        for rank, key in enumerate(sample_data.keys()):
            if key not in ranks:
                ranks[key] = []
            ranks[key].append(rank)
    ranks = [(k, mean(rs)) for k, rs in ranks.items()]
    return [key for key, _ in sorted(ranks, key=lambda x: x[1])]


def write_table(data, columns, output=None):
    with ExitStack() as estack:
        if output is None:
            ofile = sys.stdout
        else:
            ofile = estack.enter_context(open(output, 'w'))

        columns = ['sample', *columns]
        ofile.write('\t'.join(columns))
        ofile.write('\n')
        for sample_id, sample_data in data:
            sample_data['sample'] = sample_id
            row = [str(sample_data.get(col, '')) for col in columns]
            ofile.write('\t'.join(row))
            ofile.write('\n')
        if output is not None:
            print(f'[OK] written {output}')


def write_output(data, output):
    out_txt = json.dumps(data, indent=4)
    if output is None:
        print(out_txt)
    else:
        with open(output, 'w') as ofile:
            ofile.write(out_txt)
            ofile.write('\n')
        print(f'[OK] written: {output}')


def aggregate(data, dataset=None):
    GROUP_KEY = ('layout', 'model_name', 'fwd_primer', 'rev_primer')
    groups = {}
    for sample_id, samp_data in data:
        key = tuple(samp_data.pop(i) for i in GROUP_KEY)
        if key not in groups:
            groups[key] = []

        groups[key].append(dict(sample_id=sample_id, **samp_data))

    # can't json-serialize the group key tuple as mapping key (must be string?)
    json_compat_data = []
    for key, group in groups.items():
        item = dict(zip(GROUP_KEY, key))
        item['samples'] = group
        json_compat_data.append(item)
    return json_compat_data


if __name__ == '__main__':
    cli()
