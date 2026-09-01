#!/usr/bin/env python3
"""
Run MS2Query on a single MGF file (offline, no account) using pre-downloaded
Zenodo models. Run inside the ms2query conda env.
"""
import argparse
from pathlib import Path
import shutil


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

    # run_complete_folder processes every .mgf in a directory; use an
    # isolated directory containing just this one file so it can't pick up
    # unrelated MGFs from a shared parent directory.
    spectra_dir = outdir / 'ms2query_input'
    spectra_dir.mkdir(exist_ok=True)
    staged_mgf = spectra_dir / mgf.name
    if staged_mgf.exists():
        staged_mgf.unlink()
    staged_mgf.symlink_to(mgf.resolve())

    library = create_library_object_from_one_dir(args.models_dir)
    run_complete_folder(library, str(spectra_dir))

    # run_complete_folder writes to spectra_dir/results/<mgf-stem>.csv
    produced = spectra_dir / 'results' / f'{mgf.stem}.csv'
    dest = outdir / 'results.csv'
    if not produced.exists():
        raise FileNotFoundError(f'MS2Query did not produce expected output: {produced}')
    shutil.move(str(produced), str(dest))
    print(f'[OK] wrote {dest}')


if __name__ == '__main__':
    main()
