"""
Decide whether mzmine's scan-based RT correction should be applied to a
processing group, and generate the mzmine sample-metadata TSV it needs.

Real, validated finding (set_41, mzmine 4.10.6, full campaign vs. subsets of
1/2/3/4/5 acquisition batches, measured via a real drift-split diagnostic on
the aligned feature table):

    batches   uncorrected rate   corrected rate   relative change
    1         0.13%              0.15%            +15% (worse)
    2         1.44%              1.51%            +5%  (worse)
    3         1.79%              1.72%            -4%  (marginal)
    4         6.46%              4.78%            -26% (substantial)
    5         6.9%               4.2%             -39% (substantial)

RT correction interpolates between whatever calibration runs (QC-type
samples) are available *within one mzmine invocation* - with only 1-2
batches it has too few temporal knots to track real drift and net makes
things slightly worse; from 3 batches on it stops hurting and by 4+ it's a
clear, substantial win. min_batches=3 (the point it stops being harmful) is
the currently configured threshold - see config.yaml's
metabolomics.rt_correction.min_batches.

This module does NOT try to detect batches itself - that's
run_manifest.py's batch_id_guess/batch_id_override job, done once per
dataset. This module just counts distinct batch_id values within a
processing group's manifest rows and checks each one has enough QC pool
coverage to actually contribute calibration points.
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

from ..utils import logme, UsageError

DEFAULT_MIN_BATCHES = 3
DEFAULT_MIN_QC_PER_BATCH = 3

# mzmine's own metadata-file column names (see ProjectMetadataReader /
# MetadataColumn - the internal key names, distinct from the human-readable
# title row; confirmed by real testing that only these exact strings on the
# THIRD header row are recognized).
MZMINE_FILENAME_KEY = 'filename'
MZMINE_SAMPLE_TYPE_KEY = 'mzmine_sample_type'
MZMINE_RUN_DATE_KEY = 'run_date'

# ms_role -> mzmine SampleType string. Roles not listed pass through
# unchanged; SampleType.ofString() regex-matches "qc" as a substring, so
# qc_spike is correctly picked up as a QC type too without an explicit
# mapping (confirmed by real testing - see rt_correction.py's module
# docstring numbers, which already reflect this).
ROLE_TO_MZMINE_SAMPLE_TYPE = {
    'qc_pool': 'qc',
    'cal_standard': 'calibration',
}


@logme()
def eligible(manifest_rows, min_batches=DEFAULT_MIN_BATCHES,
             min_qc_per_batch=DEFAULT_MIN_QC_PER_BATCH, log=None):
    """
    manifest_rows: this processing group's rows from run_manifest.tsv (list
    of dicts, as read by csv.DictReader).

    Returns (bool, str) - the decision and a one-line human-readable reason,
    suitable for logging into the rule's log file so the choice is
    auditable after the fact.
    """
    by_batch = defaultdict(list)
    undated = 0
    for row in manifest_rows:
        batch_id = row.get('batch_id')
        if not batch_id:
            undated += 1
            continue
        by_batch[batch_id].append(row)

    n_batches = len(by_batch)
    if n_batches < min_batches:
        return False, (
            f'{n_batches} batch(es) < min_batches={min_batches} '
            f'({undated} row(s) had no batch_id) - skipping RT correction'
        )

    qc_counts = {
        batch_id: sum(1 for r in rows if r.get('ms_role') == 'qc_pool')
        for batch_id, rows in by_batch.items()
    }
    thin_batches = {b: c for b, c in qc_counts.items() if c < min_qc_per_batch}
    if thin_batches:
        return False, (
            f'{n_batches} batch(es) meet min_batches={min_batches}, but '
            f'{len(thin_batches)} have fewer than min_qc_per_batch='
            f'{min_qc_per_batch} qc_pool samples: {thin_batches} - '
            f'skipping RT correction'
        )

    return True, (
        f'{n_batches} batch(es), all with >= {min_qc_per_batch} qc_pool '
        f'samples ({qc_counts}) - applying RT correction'
    )


def write_mzmine_metadata(manifest_rows, output):
    """
    Write an mzmine project-metadata TSV for this processing group's
    samples, in the 3-header-row wide format mzmine's ProjectMetadataReader
    expects (see io.github.mzmine.modules.visualization.projectmetadata.io.
    ProjectMetadataReader.extractRemoveHeaderTitlesAndTypes).

    Row 1 is a human-readable description row - confirmed by real testing
    that its first cell must NOT be exactly "filename"/"filenames"
    (case-insensitive), or mzmine's parser mistakes it for the real title
    row below and never finds the actual column keys, silently disabling
    sample-type lookup (it falls back to regex-guessing the type from the
    raw filename instead, with no warning). "File name" (with a space, as
    used in mzmine's own norm_rtcalibration2 test fixtures) avoids this.
    """
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)

    rows_out = []
    for row in manifest_rows:
        acq = row.get('acquisition_date', '')
        if not acq:
            continue
        # manifest stores "%m/%d/%Y %H:%M:%S" (see run_manifest.py); mzmine's
        # DATETIME columns expect ISO 8601.
        from datetime import datetime
        try:
            dt = datetime.strptime(acq, '%m/%d/%Y %H:%M:%S')
        except ValueError:
            continue
        sample_type = ROLE_TO_MZMINE_SAMPLE_TYPE.get(
            row.get('ms_role', ''), row.get('ms_role', ''))
        rows_out.append((
            f"{row['sample']}.mzML",
            sample_type,
            dt.strftime('%Y-%m-%dT%H:%M:%S'),
        ))
    rows_out.sort()

    with open(output, 'w', newline='') as ofile:
        writer = csv.writer(ofile, delimiter='\t')
        writer.writerow(['File name', 'Sample type', 'Run start time stamp of the sample'])
        writer.writerow(['TEXT', 'TEXT', 'DATETIME'])
        writer.writerow([MZMINE_FILENAME_KEY, MZMINE_SAMPLE_TYPE_KEY, MZMINE_RUN_DATE_KEY])
        writer.writerows(rows_out)
    return len(rows_out)


def cli():
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument('manifest', help='run_manifest.tsv path')
    argp.add_argument('group_id', help='processing group to check/export')
    argp.add_argument('--min-batches', type=int, default=DEFAULT_MIN_BATCHES)
    argp.add_argument('--min-qc-per-batch', type=int, default=DEFAULT_MIN_QC_PER_BATCH)
    argp.add_argument('--write-metadata', metavar='OUTPUT_TSV',
                       help='also write the mzmine metadata TSV for this group')
    args = argp.parse_args()

    groups_tsv = Path(args.manifest).parent / 'processing_groups.tsv'
    with open(groups_tsv) as ifile:
        group_samples = {
            r['sample'] for r in csv.DictReader(ifile, delimiter='\t')
            if r['group_id'] == args.group_id
        }
    with open(args.manifest) as ifile:
        rows = [r for r in csv.DictReader(ifile, delimiter='\t')
                if r['sample'] in group_samples]

    ok, reason = eligible(rows, args.min_batches, args.min_qc_per_batch)
    print(f'{"ELIGIBLE" if ok else "NOT ELIGIBLE"}: {reason}')

    if args.write_metadata:
        n = write_mzmine_metadata(rows, args.write_metadata)
        print(f'wrote {n} rows to {args.write_metadata}')


if __name__ == '__main__':
    cli()
