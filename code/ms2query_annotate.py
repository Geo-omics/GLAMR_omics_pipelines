#!/usr/bin/env python3
"""
Run MS2Query on a single MGF file (offline, no account) using pre-downloaded
Zenodo models. Run inside the ms2query conda env.

The input is MZmine's SIRIUS-style MGF, which holds an MS1 (isotope pattern)
entry and one or more MS2 entries per feature. Only the MS2 entries are
queried - MS2Query would otherwise annotate the MS1 entries as if they were
fragment spectra - and a feature_id column (from FEATURE_ID=) is added to the
results so they join to the feature table without relying on MGF entry order.
query_spectrum_nr is the 1-based index into the MS2-only MGF.
"""
import argparse
import csv
from pathlib import Path
import shutil


def read_mgf_entries(path):
    """ MGF as a list of entries, each a list of lines (BEGIN/END IONS included) """
    entries, cur = [], None
    with open(path) as fh:
        for line in fh:
            if line.strip() == 'BEGIN IONS':
                cur = [line]
            elif cur is not None:
                cur.append(line)
                if line.strip() == 'END IONS':
                    entries.append(cur)
                    cur = None
    return entries


def entry_field(entry, key):
    for line in entry:
        if line.startswith(key + '='):
            return line.strip().split('=', 1)[1]
    return None


def write_ms2_mgf(mgf, dest):
    """ Writes only the MSLEVEL=2 entries of mgf to dest; returns their feature ids in order """
    ms2 = [e for e in read_mgf_entries(mgf) if entry_field(e, 'MSLEVEL') == '2']
    with open(dest, 'w') as fh:
        for e in ms2:
            fh.writelines(e)
            fh.write('\n')
    return [entry_field(e, 'FEATURE_ID') for e in ms2]


def add_feature_ids(results_csv, feature_ids, dest):
    """ Adds feature_id (looked up by 1-based query_spectrum_nr) as the first column """
    with open(results_csv, newline='') as fh:
        rows = list(csv.DictReader(fh))
    fields = ['feature_id'] + (list(rows[0].keys()) if rows else [])
    with open(dest, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({'feature_id': feature_ids[int(r['query_spectrum_nr']) - 1], **r})


def main():
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument('--mgf', required=True)
    argp.add_argument('--models-dir', required=True)
    argp.add_argument('--outdir', required=True)
    args = argp.parse_args()

    from ms2query.ms2library import create_library_object_from_one_dir
    from ms2query.run_ms2query import run_complete_folder

    mgf = Path(args.mgf)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # run_complete_folder processes every .mgf in a directory, and if its
    # results file already exists (e.g. left by a killed run) it writes to
    # "<stem>(1).csv" instead - so start from an empty, isolated directory.
    spectra_dir = outdir / 'ms2query_input'
    shutil.rmtree(spectra_dir, ignore_errors=True)
    spectra_dir.mkdir()
    staged_mgf = spectra_dir / f'{mgf.stem}_ms2.mgf'
    feature_ids = write_ms2_mgf(mgf, staged_mgf)
    print(f'[INFO] {len(feature_ids)} MS2 spectra to query')

    library = create_library_object_from_one_dir(args.models_dir)
    run_complete_folder(library, str(spectra_dir))

    produced = spectra_dir / 'results' / f'{staged_mgf.stem}.csv'
    if not produced.exists():
        raise FileNotFoundError(f'MS2Query did not produce expected output: {produced}')
    dest = outdir / 'results.csv'
    add_feature_ids(produced, feature_ids, dest)
    staged_mgf.unlink()
    print(f'[OK] wrote {dest}')


if __name__ == '__main__':
    main()
