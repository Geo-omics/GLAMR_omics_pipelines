#!/usr/bin/env python3
"""
Download MS2Query's pretrained library/model files from Zenodo, once, for a
given ionization mode. Run inside the ms2query conda env.
"""
import argparse
from pathlib import Path


def main():
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument('--ion-mode', required=True, choices=['positive', 'negative'])
    argp.add_argument('--outdir', required=True)
    argp.add_argument('--marker', required=True, help='marker file to touch on success')
    args = argp.parse_args()

    from ms2query.run_ms2query import download_zenodo_files

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    download_zenodo_files(args.ion_mode, str(outdir))

    Path(args.marker).touch()
    print(f'[OK] downloaded MS2Query {args.ion_mode} models to {outdir}')


if __name__ == '__main__':
    main()
