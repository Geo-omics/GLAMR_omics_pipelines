"""
Assign metabolomics samples to processing groups.

A processing group is the unit of MZmine alignment: (method_profile x
polarity). Alignment across dissimilar methods or polarities is invalid, so
grouping by these two fields - not just by dataset - is what lets a single
dataset that mixes e.g. positive and negative mode, or two different
instruments, still produce correct per-group feature matrices. For a
homogeneous project like set_41 this collapses to exactly one group.
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

from ..utils import logme, UsageError


def cli():
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument('manifest', help='run_manifest.tsv path')
    argp.add_argument('method_profile', help="the dataset's assigned method_profile name")
    argp.add_argument('--output', default='processing_groups.tsv')
    args = argp.parse_args()
    try:
        with open(args.manifest) as ifile:
            rows = list(csv.DictReader(ifile, delimiter='\t'))
        groups = assign(rows, args.method_profile)
        write_groups(groups, args.output)
    except UsageError as e:
        argp.error(e)


@logme()
def assign(manifest_rows, method_profile, log=None):
    """
    Partition manifest rows into processing groups.

    manifest_rows: list of dicts as read from run_manifest.tsv
    method_profile: the project's assigned method_profile name (one dataset
        is assumed to use one instrument method; a dataset genuinely
        spanning multiple methods should be split into separate StudyIDs
        upstream, matching how GLAMR already scopes StudyID to one
        acquisition campaign)

    Returns: dict of group_id -> list of sample manifest rows. All roles
    (sample/blank/qc_pool/qc_spike/cal_standard/wash) are grouped together -
    QC-based filtering happens downstream (stage 5), not here.
    """
    groups = defaultdict(list)
    for row in manifest_rows:
        polarity = row.get('polarity') or ''
        if polarity not in ('positive', 'negative'):
            # Polarity-switching runs must be split into per-polarity mzML
            # at the conversion stage (stage 1) before they ever reach
            # manifest/group assignment. If one slipped through, fail loudly
            # rather than silently building a mixed-polarity group, which
            # would silently corrupt alignment.
            raise UsageError(
                f'sample {row["sample"]!r} has no single polarity '
                f'(got {polarity!r}) - polarity-switching runs must be '
                f'split into per-polarity mzML before group assignment'
            )
        group_id = f'{method_profile}_{polarity}'
        groups[group_id].append(row)
    return dict(groups)


def write_groups(groups, output):
    """ Write group_id, polarity, n_samples, sample (one row per sample) """
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, 'w', newline='') as ofile:
        writer = csv.writer(ofile, delimiter='\t')
        writer.writerow(['group_id', 'sample', 'ms_role'])
        for group_id, rows in sorted(groups.items()):
            for row in sorted(rows, key=lambda r: _sample_num(r['sample'])):
                writer.writerow([group_id, row['sample'], row['ms_role']])
    print(f'[OK] wrote {output}: {len(groups)} group(s), '
          f'{sum(len(r) for r in groups.values())} samples')


def _sample_num(sample):
    try:
        return int(sample.removeprefix('samp_'))
    except ValueError:
        return sample


def get_group_samples(groups_tsv, group_id):
    """ Input-function helper: list of sample IDs belonging to a group """
    samples = []
    with open(groups_tsv) as ifile:
        for row in csv.DictReader(ifile, delimiter='\t'):
            if row['group_id'] == group_id:
                samples.append(row['sample'])
    return samples


def get_all_group_ids(groups_tsv):
    """ Input-function helper: list of all group IDs in a groups file """
    ids = []
    with open(groups_tsv) as ifile:
        for row in csv.DictReader(ifile, delimiter='\t'):
            if row['group_id'] not in ids:
                ids.append(row['group_id'])
    return ids


if __name__ == '__main__':
    cli()
