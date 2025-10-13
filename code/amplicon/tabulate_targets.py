"""
Collect targeting information for all samples of a dataset
"""
import argparse
from collections import Counter
from pathlib import Path
import re

from . import Err


TEST_FILE_PATTERN = 'guess*{sample_id}.txt'


def main():
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
        '--test-dir',
        help='Directory with target info files for testing.  Target info files'
             f'must be named via this pattern: "{TEST_FILE_PATTERN}"',
    )
    args = argp.parse_args()

    project_dir = Path(args.project_directory)

    data = collect_target_data(project_dir, args.test_dir)
    data = mangle_data(data)
    print_preamble(data, project_dir.name)
    print_table(data)


def collect_target_data(project_dir, test_dir=None):
    """
    Collect data from target_info files under project directory
    """
    data = []
    for i in (project_dir / 'amplicons').iterdir():
        if not i.is_dir():
            continue

        if i.is_symlink() and i.name != i.resolve().name:
            # TODO / FIXME - meant for dev/testing only
            # interpret this as skip / do not process sample
            continue
        sample_id = i.name
        if test_dir:
            test_file_glob = TEST_FILE_PATTERN.format(sample_id=sample_id)
            for j in Path(test_dir).glob(test_file_glob):
                input_file = j
                break
            else:
                input_file = None
        else:
            input_file = i / 'detect_region' / 'target_info.txt'

        if input_file:
            row = get_target_info(input_file)
        else:
            row = {'errors': [(Err.E19, 'no test target info file')]}
        row['sample_id'] = sample_id
        data.append(row)

    return data


error_pat = re.compile(r'^error E([0-9]+)')


def get_target_info(target_info_file):
    """
    Read target info file for a sample
    """
    data = {'errors': []}
    err_map = {i.value: i for i in Err}
    try:
        ifile = open(target_info_file)
    except Exception as e:
        data['errors'].append(
            (Err.E18, f'no target info: {e.__class__.__name__}: {e}')
        )
        return data

    with ifile as ifile:
        for lnum, line in enumerate(ifile, start=1):
            line = line.strip()
            key, _, value = line.partition(':')
            value = value.strip()
            try:
                value = int(value)
            except ValueError:
                pass

            if key in data:
                raise ValueError(
                    f'duplicate key in {target_info_file} line {lnum}'
                )
            if m := error_pat.match(key):
                try:
                    err = err_map[int(m.group(1))]
                except KeyError as e:
                    raise ValueError('invalid error number') from e
                data['errors'].append((err, value))
            else:
                data[key] = value
    return data


def mangle_data(data):
    for row in data:
        row['errors_full'] = row['errors']
        row['errors'] = ','.join(err.name for err, _ in row['errors_full'])
        if row.get('swapped') == 'YES':
            row['swapped'] = 'swapped'
    return data


def print_preamble(data, dataset):
    errors = Counter()
    for row in data:
        errors.update(row['errors_full'])

    print(f'# dataset {dataset} with {len(data)} amplicon samples')
    cats = Counter()
    for row in data:
        if row['errors']:
            cats['errors'] += 1
        elif row['primers']:
            cats[row['primers']] += 1
        else:
            raise ValueError('expecting primer or error')
    print('#', *cats.most_common())

    if errors:
        for (err, msg), count in sorted(errors.items(), key=lambda x: x[0][0]):
            print(f'# error {err.name}: [{count}x]', msg)
    else:
        print('# no errors')


def print_table(data):
    columns = (
        'sample_id', 'errors', 'swapped', 'primers', 'model', 'primer_content',
        'distance', 'length', 'overlap',
    )

    def sortkey(row):
        try:
            num = int(row['sample_id'].removeprefix('samp_'))
        except Exception:
            num = 0
        return (row.get('model', ''), row.get('primers', ''), num)

    print(*columns, sep='\t')
    for row in sorted(data, key=sortkey):
        print(*(row.get(col, '') for col in columns), sep='\t')


if __name__ == '__main__':
    main()
