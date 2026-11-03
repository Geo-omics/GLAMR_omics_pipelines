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

(Those numbers came from an earlier diagnostic whose code was not kept; they
are not comparable with pypelib.metabolomics.drift_split, which is the
diagnostic to use from now on.)

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

Reference runs ("qc" mode vs "samples" mode)
--------------------------------------------
mzmine's ScanRtCorrectionTask only uses landmark features ("standards") that
are found in EVERY reference run (RtStandard.isValid() - confirmed by reading
the v4.10.6 source), and corrects all other runs by interpolating between the
nearest reference runs in acquisition time. In "qc" mode the references are
the qc_pool/qc_spike runs. On set_41 that under-corrected badly: the QC pools
share few compounds with the real samples, and the measured between-batch
shift is strongly RT-dependent (raw apex RTs, 2024-06-20 vs 06-25 batches:
~0 min at 2-5 min, +0.28 min at 7-8 min, ~0 again at 12-13 min). MC-LR still
ended up split across five aligned features after "qc"-mode correction.

"samples" mode instead picks a few real sample runs per batch as references
(one per injection-order bin, preferring runs with the most MS2 scans so that
failed injections are avoided), labels them mzmine sample type "qc", and
labels the real QC runs "sample" so they no longer constrain the landmark set.
The mzmine sample type is only used by this RT-correction step - downstream
QC filtering uses the manifest's ms_role, which is unchanged.
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

from ..utils import logme, UsageError

DEFAULT_MIN_BATCHES = 3
DEFAULT_MIN_QC_PER_BATCH = 3
DEFAULT_REFERENCE_MODE = 'qc'
DEFAULT_REFERENCES_PER_BATCH = 4
REFERENCE_MODES = ('qc', 'samples')

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


def _reference_role(reference_mode):
    """ The ms_role whose runs can serve as RT references in this mode """
    if reference_mode not in REFERENCE_MODES:
        raise UsageError(f'unknown RT reference mode {reference_mode!r}, expected one of {REFERENCE_MODES}')
    return 'qc_pool' if reference_mode == 'qc' else 'sample'


def select_reference_samples(manifest_rows, per_batch=DEFAULT_REFERENCES_PER_BATCH, quality=None):
    """
    Pick up to per_batch real sample runs per batch to act as RT-correction
    references: split each batch's sample runs into per_batch equal
    injection-order bins and take the highest-quality run from each bin.

    quality: optional dict sample -> number (higher is better; the pipeline
    passes the MS2 scan count from run_metadata.json, which is lowest for
    failed injections). Without it, the middle run of each bin is used.

    Returns the set of selected sample IDs.
    """
    from datetime import datetime
    quality = quality or {}
    by_batch = defaultdict(list)
    for row in manifest_rows:
        if row.get('ms_role') != 'sample' or not row.get('batch_id'):
            continue
        try:
            acq = datetime.strptime(row.get('acquisition_date', ''), '%m/%d/%Y %H:%M:%S')
        except ValueError:
            continue
        by_batch[row['batch_id']].append((acq, row['sample']))

    selected = set()
    for runs in by_batch.values():
        runs.sort()
        n_bins = min(per_batch, len(runs))
        for b in range(n_bins):
            bin_runs = runs[b * len(runs) // n_bins:(b + 1) * len(runs) // n_bins]
            if quality:
                best = max(bin_runs, key=lambda r: quality.get(r[1], float('-inf')))
            else:
                best = bin_runs[len(bin_runs) // 2]
            selected.add(best[1])
    return selected


@logme()
def eligible(manifest_rows, min_batches=DEFAULT_MIN_BATCHES,
             min_qc_per_batch=DEFAULT_MIN_QC_PER_BATCH,
             reference_mode=DEFAULT_REFERENCE_MODE, log=None):
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

    # min_qc_per_batch applies to whichever role supplies the references
    ref_role = _reference_role(reference_mode)
    qc_counts = {
        batch_id: sum(1 for r in rows if r.get('ms_role') == ref_role)
        for batch_id, rows in by_batch.items()
    }
    thin_batches = {b: c for b, c in qc_counts.items() if c < min_qc_per_batch}
    if thin_batches:
        return False, (
            f'{n_batches} batch(es) meet min_batches={min_batches}, but '
            f'{len(thin_batches)} have fewer than min_qc_per_batch='
            f'{min_qc_per_batch} {ref_role} runs: {thin_batches} - '
            f'skipping RT correction'
        )

    return True, (
        f'{n_batches} batch(es), all with >= {min_qc_per_batch} {ref_role} '
        f'runs ({qc_counts}) - applying RT correction ({reference_mode!r} references)'
    )


def write_mzmine_metadata(manifest_rows, output, reference_mode=DEFAULT_REFERENCE_MODE,
                          references_per_batch=DEFAULT_REFERENCES_PER_BATCH, quality=None):
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

    reference_mode: "qc" labels qc_pool runs as mzmine "qc" (the RT-correction
    references); "samples" labels select_reference_samples() runs as "qc"
    and every other non-blank, non-calibration run as "sample" - see the
    module docstring. Returns (rows written, set of reference sample IDs).
    """
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)

    references = (select_reference_samples(manifest_rows, references_per_batch, quality)
                  if _reference_role(reference_mode) == 'sample' else None)

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
        role = row.get('ms_role', '')
        if references is None:
            sample_type = ROLE_TO_MZMINE_SAMPLE_TYPE.get(role, role)
        elif row['sample'] in references:
            sample_type = 'qc'
        elif role in ('blank', 'cal_standard'):
            sample_type = ROLE_TO_MZMINE_SAMPLE_TYPE.get(role, role)
        else:
            # includes qc_pool / qc_spike: any string containing "qc" would make
            # mzmine treat the run as a reference (SampleType.ofString)
            sample_type = 'sample'
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
    if references is None:
        references = {r['sample'] for r in manifest_rows if r.get('ms_role') == 'qc_pool'}
    return len(rows_out), references


def cli():
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument('manifest', help='run_manifest.tsv path')
    argp.add_argument('group_id', help='processing group to check/export')
    argp.add_argument('--min-batches', type=int, default=DEFAULT_MIN_BATCHES)
    argp.add_argument('--min-qc-per-batch', type=int, default=DEFAULT_MIN_QC_PER_BATCH)
    argp.add_argument('--reference-mode', choices=REFERENCE_MODES, default=DEFAULT_REFERENCE_MODE)
    argp.add_argument('--references-per-batch', type=int, default=DEFAULT_REFERENCES_PER_BATCH)
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

    ok, reason = eligible(rows, args.min_batches, args.min_qc_per_batch, args.reference_mode)
    print(f'{"ELIGIBLE" if ok else "NOT ELIGIBLE"}: {reason}')

    if args.write_metadata:
        n, refs = write_mzmine_metadata(rows, args.write_metadata, args.reference_mode,
                                        args.references_per_batch)
        print(f'wrote {n} rows ({len(refs)} RT references) to {args.write_metadata}')


if __name__ == '__main__':
    cli()
