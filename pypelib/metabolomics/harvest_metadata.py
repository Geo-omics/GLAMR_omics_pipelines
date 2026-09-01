"""
Harvest instrument/acquisition metadata from a converted mzML file, with
optional enrichment from ThermoRawFileParser's metadata-only JSON mode when
the source was a Thermo raw file.

Everything here is auto-derived. The only fields the pipeline cannot derive
and must be told by a human are recorded separately in the run manifest
(see run_manifest.py): ms_role, prep_variant, cal_level, BioSampleID,
amount_basis.

NOTE: ThermoRawFileParser's metadata-only mode (-f 4 -m 0) does NOT expose
scan polarity - confirmed by testing against real GLAMR data (its JSON has
FileProperties/InstrumentProperties/MsData/ScanSettings/SampleData sections,
none of which carry polarity). So unlike what a metadata-only fast path
would suggest, polarity/MS-levels/centroid-vs-profile must come from the
converted mzML - this rule runs after conversion for every vendor, Thermo
included. The Thermo JSON is still useful for richer fields (instrument
model, method file path, sample name) and is merged in when available.
"""
import argparse
import json
from pathlib import Path
from xml.etree import ElementTree as ET

from ..utils import conda_run, logme, UsageError

MZML_NS = '{http://psi.hupo.org/ms/mzml}'


def cli():
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument('mzml', help='Path to the converted mzML file')
    argp.add_argument('output', help='Path to write the metadata JSON to')
    argp.add_argument(
        '--thermo-raw',
        help='Path to the original Thermo .raw file, for enrichment (optional)',
    )
    argp.add_argument('--conda-env', help='conda env prefix with ThermoRawFileParser')
    args = argp.parse_args()
    try:
        main(args.mzml, args.output, thermo_raw=args.thermo_raw, conda_env=args.conda_env)
    except UsageError as e:
        argp.error(e)


@logme()
def main(mzml_path, output, thermo_raw=None, conda_env=None, log=None):
    metadata = _harvest_mzml(Path(mzml_path))

    if thermo_raw is not None:
        try:
            enrichment = _harvest_thermo(Path(thermo_raw), conda_env)
        except Exception as e:
            print(f'[WARNING] Thermo metadata enrichment failed, continuing '
                  f'without it: {e.__class__.__name__}: {e}')
        else:
            # mzML-derived fields (polarity, ms levels, format) stay
            # authoritative on conflicts; missing/None mzML values get
            # filled in from the Thermo enrichment.
            for k, v in enrichment.items():
                if metadata.get(k) is None:
                    metadata[k] = v

    metadata['source_mzml'] = str(mzml_path)

    with open(output, 'w') as ofile:
        json.dump(metadata, ofile, indent=4, sort_keys=True)
        ofile.write('\n')
    print(f'[OK] wrote {output}')
    return metadata


def _harvest_mzml(mzml_path):
    """
    Pull polarity / MS levels / centroid-vs-profile / scan count directly out
    of the mzML header via streaming XML parse, no extra dependency needed.
    """
    polarities = set()
    ms_levels = set()
    formats = set()
    scan_count = 0
    instrument_model = None

    for event, elem in ET.iterparse(str(mzml_path), events=('end',)):
        tag = elem.tag.removeprefix(MZML_NS)
        if tag == 'spectrum':
            scan_count += 1
            for cv in elem.iter(f'{MZML_NS}cvParam'):
                name = cv.get('name', '')
                if name == 'positive scan':
                    polarities.add('positive')
                elif name == 'negative scan':
                    polarities.add('negative')
                elif name == 'ms level':
                    ms_levels.add(int(cv.get('value')))
                elif name == 'centroid spectrum':
                    formats.add('centroid')
                elif name == 'profile spectrum':
                    formats.add('profile')
            elem.clear()
        elif tag in ('referenceableParamGroup', 'instrumentConfiguration'):
            for cv in elem.iter(f'{MZML_NS}cvParam'):
                if instrument_model is None and 'instrument model' in cv.get('name', '').lower():
                    instrument_model = cv.get('name')

    if scan_count == 0:
        raise UsageError(f'no spectra found in {mzml_path} - conversion may have failed')

    if len(polarities) > 1:
        acquisition_mode = 'polarity_switching'
    elif len(polarities) == 1:
        acquisition_mode = None  # single polarity, nothing special to flag
    else:
        acquisition_mode = None

    return {
        'instrument_model': instrument_model,
        'polarity': sorted(polarities) if len(polarities) != 1 else next(iter(polarities)),
        'polarity_switching': len(polarities) > 1,
        'ms_levels_present': sorted(ms_levels),
        'ms_format': sorted(formats) if len(formats) != 1 else (next(iter(formats)) if formats else None),
        'scan_count': scan_count,
    }


def _harvest_thermo(raw_path, conda_env):
    """
    Enrichment only (not authoritative for polarity/ms-levels/format - see
    module docstring). Uses ThermoRawFileParser's metadata-only mode:
    -f 4 skips spectra output entirely (fast, no mzML written), -m 0 writes
    JSON metadata to the explicit path given by -c.
    """
    import tempfile
    with tempfile.TemporaryDirectory() as tmpd:
        out_json = Path(tmpd) / 'meta.json'
        cmd = [
            'ThermoRawFileParser',
            '-i', str(raw_path),
            '-f', '4',
            '-m', '0',
            '-c', str(out_json),
        ]
        conda_run(conda_env, cmd, check=True)
        if not out_json.exists():
            raise UsageError(f'ThermoRawFileParser produced no metadata file for {raw_path}')
        with open(out_json) as ifile:
            trfp_meta = json.load(ifile)

    kv = {}
    for rows in trfp_meta.values():
        for row in rows:
            kv[row.get('name')] = row.get('value')

    return {
        'instrument_model': kv.get('Thermo Scientific instrument model'),
        'instrument_method': kv.get('device acquisition method'),
        'sample_name': kv.get('sample name'),
        'acquisition_date': kv.get('Content Creation Date'),
        'thermo_scan_count': kv.get('Number of scans'),
        'thermo_ms1_scans': kv.get('Number of MS1 spectra'),
        'thermo_ms2_scans': kv.get('Number of MS2 spectra'),
        'thermo_rt_range_min': kv.get('Retention time range'),
        'thermo_mz_range': kv.get('Mz range'),
        'fragmentation_type': kv.get('beam-type collision-induced dissociation') and 'HCD',
    }


if __name__ == '__main__':
    cli()
