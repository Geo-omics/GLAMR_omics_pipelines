#!/usr/bin/env python3
"""
CyanoMetDB as input for a SIRIUS custom structure database (rule
metab_sirius_cyanometdb). Reads the CyanoMetDB spreadsheet bundled with the
CPM package (the same version CPM matches against) and writes:

  --smiles   SIRIUS `custom-db import` input: tab-separated SMILES, id, name,
             no header
  --lookup   one row per compound (id, name, class, formula, masses, InChIKey,
             NPAtlas id, producer) for joining SIRIUS hits back to CyanoMetDB.
             SIRIUS reports a custom-database hit under the structure's PubChem
             name and links, so join on InChIkey2D (first InChIKey block).

Run inside the cpm conda env (pandas, openpyxl, cpm).
"""
import argparse
import importlib.resources

import pandas as pd

SHEET = 'CyanoMetDB_V03_2024'
XLSX = 'data/CyanoMetDB_Version03.xlsx'


def clean(s):
    """ Single-line text, safe in a TSV """
    return s.fillna('').astype(str).str.replace(r'[\t\r\n]+', ' ', regex=True).str.strip()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--smiles', required=True)
    ap.add_argument('--lookup', required=True)
    args = ap.parse_args()

    xlsx = importlib.resources.files('cpm') / XLSX
    db = pd.read_excel(xlsx, sheet_name=SHEET, dtype=str, engine='openpyxl')
    db.columns = db.columns.str.strip()  # v03 has 'Compound identifier ' (trailing space)
    print(f'[INFO] {len(db)} compounds in {xlsx} ({SHEET})')

    out = pd.DataFrame({
        'id': clean(db['Compound identifier']),
        'name': clean(db['Compound name']),
        'class': clean(db['Class of compound']),
        'formula': clean(db['Molecular formula (neutral)']),
        'monoisotopic_mass': clean(db['Monoisotopic mass (neutral)']),
        'mz_mh': clean(db['Monoisotopic mass [M+H]+']),
        'mz_mmh': clean(db['Monoisotopic mass [M-H]-']),
        'smiles': clean(db['SMILES (canonical or isomeric)']),
        'inchikey': clean(db['InChlKey']),  # sic, column name in v03
        'npatlas_id': clean(db['NPAtlas_ID https://www.npatlas.org']),
        'genus': clean(db['Genus']),
    })
    out['inchikey2d'] = out['inchikey'].str.split('-').str[0]
    out = out[out['id'] != '']
    has_smiles = out['smiles'] != ''
    print(f'[INFO] {has_smiles.sum()} with SMILES, {(~has_smiles).sum()} without (skipped for SIRIUS)')

    out.to_csv(args.lookup, sep='\t', index=False)
    out.loc[has_smiles, ['smiles', 'id', 'name']].to_csv(args.smiles, sep='\t', index=False, header=False)


if __name__ == '__main__':
    main()
