"""
Build the per-dataset metabolomics run manifest.

Combines auto-harvested per-sample metadata (data/omics/metabolomes/{sample}/
spectra/run_metadata.json, written by harvest_metadata.py) with an
auto-guessed run role (sample/blank/qc_pool/qc_spike/cal_standard/wash), an
auto-guessed acquisition batch_id, and a handful of fields that genuinely
cannot be derived and must be supplied by a human: prep_variant, cal_level,
BioSampleID, amount_basis.

Mirrors pypelib.amplicon.dispatch's guess/override split: ms_role_guess and
batch_id_guess are regenerated every run, *_override are human-editable
columns that are always preserved verbatim across reruns and win when set.
The other human fields (prep_variant, cal_level, BioSampleID, amount_basis)
have no meaningful auto-guess at all, so they are simply carried forward
unchanged from any existing manifest.

batch_id exists for one reason: it's the unit RT-correction eligibility is
decided on (see pypelib.metabolomics.rt_correction) - real drift correction
needs multiple, well-separated acquisition sessions to have any calibration
knots to interpolate between (confirmed empirically: 1-2 batches makes
correction slightly worse, 3 is marginal, 4+ is a real win - see
rt_correction.py's docstring for the numbers). It is NOT the same thing as a
"collection year" or any other project-specific folder convention - those
happened to line up with real acquisition sessions for set_41, but won't for
every project, so batch_id is auto-detected purely from acquisition_date
gaps, independent of file paths or naming.
"""
import argparse
import csv
import filecmp
from datetime import datetime, timedelta
import json
from pathlib import Path
import re
import shutil
import statistics
from tempfile import NamedTemporaryFile

from ..utils import logme, UsageError

# Filename-pattern first guesses for run role. These match the naming
# convention observed in the set_41 (USGS/GLERL) project; other projects
# will have different conventions and are expected to need ms_role_override
# curation - this is a starting point, not ground truth.
ROLE_PATTERNS = [
    (re.compile(r'^blank[_\s-]?\d*$', re.IGNORECASE), 'blank'),
    (re.compile(r'^(wash|blk)[_\s-]?\d*$', re.IGNORECASE), 'blank'),
    (re.compile(r'^QC[_-]?SA[_-]?\d*$', re.IGNORECASE), 'qc_spike'),
    (re.compile(r'^QC[_-]?A?[_-]?\d*$', re.IGNORECASE), 'qc_pool'),
    (re.compile(r'^pool[_\s-]?\d*$', re.IGNORECASE), 'qc_pool'),
    (re.compile(r'^(std|cal|calib)', re.IGNORECASE), 'cal_standard'),
    (re.compile(r'^LCTXM', re.IGNORECASE), 'cal_standard'),
]

VALID_ROLES = {'sample', 'blank', 'qc_pool', 'qc_spike', 'cal_standard', 'wash'}

MANIFEST_COLUMNS = [
    'sample', 'ms_role', 'ms_role_guess', 'ms_role_override',
    'batch_id', 'batch_id_guess', 'batch_id_override',
    'prep_variant', 'cal_level', 'BioSampleID', 'amount_basis',
    'polarity', 'instrument_model', 'ms_levels_present', 'ms_format',
    'acquisition_date', 'source_filename',
]

# Minimum gap-detection threshold, applied even when the run's own injections
# are so close together that 10x their median gap would be implausibly small
# (e.g. a handful of samples run back-to-back). Real acquisition-session
# breaks in the data we've seen are hours to days; this floor just guards
# against pathologically small batches on tiny test datasets.
MIN_BATCH_GAP = timedelta(minutes=30)
# Multiple of the median inter-injection gap that counts as a batch boundary.
# Chosen from real data: typical within-batch injection gaps were ~15-20 min
# for set_41's method, while every real batch boundary was >=5 hours (16hr to
# ~3 days) - a gap an order of magnitude past the norm is unambiguously a new
# acquisition session, not injection-to-injection variation.
BATCH_GAP_MULTIPLIER = 10

DEFAULT_MANIFEST_NAME = 'run_manifest.tsv'


def cli():
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument(
        'project_directory',
        help='data/projects/<dataset> - expected to contain a "metabolomes" '
             'subdirectory with one directory (or symlink) per sample.',
    )
    argp.add_argument('--output', default=DEFAULT_MANIFEST_NAME)
    args = argp.parse_args()
    try:
        main(None, args.project_directory, output=args.output)
    except UsageError as e:
        argp.error(e)


@logme()
def main(metadata_files, project_dir, output=None, log=None):
    """
    metadata_files: list of run_metadata.json paths (passed by Snakemake as
    the rule's `input`, but we re-derive the sample list from the project
    directory so this also works when called directly).
    """
    project_dir = Path(project_dir).resolve()
    samples_dir = project_dir / 'metabolomes'
    if not samples_dir.is_dir():
        raise UsageError(f'no such directory: {samples_dir}')

    output = Path(output) if output else project_dir / DEFAULT_MANIFEST_NAME

    existing = _read_existing(output)

    rows = []
    for sample_dir in sorted(samples_dir.glob('samp_*'), key=_sample_sort_key):
        sample = sample_dir.name
        meta_path = sample_dir.resolve() / 'spectra' / 'run_metadata.json'
        if not meta_path.is_file():
            print(f'[WARNING] no harvested metadata for {sample}, skipping '
                  f'(expected {meta_path})')
            continue

        with open(meta_path) as ifile:
            meta = json.load(ifile)

        source_filename = _get_original_filename(sample_dir)
        guess = _guess_role(source_filename)

        prior = existing.get(sample, {})
        override = prior.get('ms_role_override', '')
        role = override if override else guess

        batch_override = prior.get('batch_id_override', '')

        rows.append({
            'sample': sample,
            'ms_role': role,
            'ms_role_guess': guess,
            'ms_role_override': override,
            'batch_id_override': batch_override,
            'prep_variant': prior.get('prep_variant', ''),
            'cal_level': prior.get('cal_level', ''),
            'BioSampleID': prior.get('BioSampleID', ''),
            'amount_basis': prior.get('amount_basis', ''),
            'polarity': meta.get('polarity', ''),
            'instrument_model': meta.get('instrument_model', ''),
            'ms_levels_present': ','.join(str(i) for i in meta.get('ms_levels_present', [])),
            'ms_format': meta.get('ms_format', ''),
            'acquisition_date': meta.get('acquisition_date', ''),
            'source_filename': source_filename,
        })

    if not rows:
        raise UsageError(f'no samples with harvested metadata found under {samples_dir}')

    _guess_batches(rows)
    for row in rows:
        row['batch_id'] = row['batch_id_override'] or row['batch_id_guess']

    _report_role_counts(rows)
    _report_batch_counts(rows)
    _write_manifest(rows, output)
    return rows


def _sample_sort_key(path):
    try:
        return int(path.name.removeprefix('samp_'))
    except ValueError:
        return path.name


def _get_original_filename(sample_dir):
    """
    Recover the original vendor filename via the import symlink, so role
    guessing works without depending on any project-specific import table.
    Works for both file-based (Thermo .raw) and directory-based
    (Agilent/Bruker .d, Waters .raw) vendor formats.
    """
    candidates = list((sample_dir / 'spectra').glob('raw_spectra.*'))
    if not candidates:
        raise UsageError(f'no raw_spectra.* found under {sample_dir / "spectra"}')
    if len(candidates) > 1:
        raise UsageError(
            f'expected exactly one raw_spectra.* under {sample_dir / "spectra"}, '
            f'found {len(candidates)}: {candidates}'
        )
    return candidates[0].resolve().name


def _guess_role(filename):
    stem = Path(filename).stem
    for pattern, role in ROLE_PATTERNS:
        if pattern.match(stem):
            return role
    return 'sample'


def _guess_batches(rows):
    """
    Auto-detect acquisition batches from real gaps in acquisition_date,
    mutating each row in place to add 'batch_id_guess'.

    Deliberately independent of any folder/filename convention: it only
    looks at acquisition timestamps, sorted dataset-wide, and calls a new
    batch wherever the gap since the previous injection is far outside the
    typical injection-to-injection interval (see BATCH_GAP_MULTIPLIER's
    docstring for why 10x the median gap was chosen). Each batch is labeled
    by the date of its first injection (e.g. "2024-06-16"), which reads
    naturally in the manifest and is stable as long as two real batches
    don't start on the same calendar day.

    Rows with no acquisition_date (harvest failed, or a vendor that doesn't
    expose one) get batch_id_guess = '' - they simply don't participate in
    batch-count-based RT-correction eligibility (see rt_correction.py).
    """
    dated = []
    for row in rows:
        row['batch_id_guess'] = ''
        raw = row.get('acquisition_date')
        if not raw:
            continue
        try:
            # real format written by harvest_metadata.py, from mzML's
            # Content Creation Date CV term (e.g. "06/17/2024 00:04:56")
            dt = datetime.strptime(raw, '%m/%d/%Y %H:%M:%S')
        except ValueError:
            continue
        dated.append((dt, row))

    dated.sort(key=lambda pair: pair[0])

    if len(dated) < 2:
        # Nothing to compare gaps against - a single dated sample is its own batch.
        for dt, row in dated:
            row['batch_id_guess'] = dt.strftime('%Y-%m-%d')
        return

    gaps = [
        (dated[i][0] - dated[i - 1][0]) for i in range(1, len(dated))
    ]
    threshold = max(MIN_BATCH_GAP, statistics.median(gaps) * BATCH_GAP_MULTIPLIER)

    batch_start = dated[0][0]
    dated[0][1]['batch_id_guess'] = batch_start.strftime('%Y-%m-%d')
    for i in range(1, len(dated)):
        if gaps[i - 1] > threshold:
            batch_start = dated[i][0]
        dated[i][1]['batch_id_guess'] = batch_start.strftime('%Y-%m-%d')


def _report_role_counts(rows):
    from collections import Counter
    counts = Counter(r['ms_role'] for r in rows)
    print(f'Role assignment for {len(rows)} samples:')
    for role, count in counts.most_common():
        print(f'  {count:>4}  {role}')
    bad_roles = {r['ms_role'] for r in rows} - VALID_ROLES
    if bad_roles:
        print(f'[WARNING] manifest has non-standard ms_role values from '
              f'manual override: {sorted(bad_roles)}')


def _report_batch_counts(rows):
    from collections import Counter
    counts = Counter(r['batch_id'] for r in rows if r['batch_id'])
    n_undated = sum(1 for r in rows if not r['batch_id'])
    print(f'Batch assignment: {len(counts)} batch(es) detected')
    for batch_id, count in sorted(counts.items()):
        print(f'  {count:>4}  {batch_id}')
    if n_undated:
        print(f'  {n_undated:>4}  (no acquisition_date - not batch-assigned)')


def _read_existing(path):
    """ Read a prior manifest, keyed by sample, for override preservation """
    if not path.is_file():
        return {}
    with open(path) as ifile:
        reader = csv.DictReader(ifile, delimiter='\t')
        return {row['sample']: row for row in reader}


def _write_manifest(rows, output_file):
    """ Idempotent write: only touches the file if content actually changed """
    output_file = Path(output_file)
    print(f'Writing {output_file} ... ', end='', flush=True)
    with NamedTemporaryFile(mode='w+', newline='') as tmpfile:
        writer = csv.DictWriter(tmpfile, fieldnames=MANIFEST_COLUMNS, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
        tmpfile.flush()
        tmpfile.seek(0)
        if output_file.exists():
            if filecmp.cmp(output_file, tmpfile.name, shallow=False):
                print('already up-to-date [OK]')
                return
            mtime = datetime.fromtimestamp(output_file.stat().st_mtime).isoformat()
            backup = output_file.with_name(output_file.name + '.' + mtime)
            print(f'making backup: {backup} ... ', end='', flush=True)
            output_file.rename(backup)
        else:
            output_file.parent.mkdir(parents=True, exist_ok=True)
        with open(output_file, 'w') as ofile:
            tmpfile.seek(0)
            shutil.copyfileobj(tmpfile, ofile)
    print('[Done]')


if __name__ == '__main__':
    cli()
