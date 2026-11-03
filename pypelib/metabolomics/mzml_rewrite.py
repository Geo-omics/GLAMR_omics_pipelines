"""
Line-by-line rewriting of indexedmzML files with the byte-offset index and
SHA-1 file checksum rebuilt, so readers that seek by index (MZmine,
pyteomics) still work after an edit.

Used by rt_align (corrected scan start times) and by rule metab_mzml_fixup,
which repairs msconvert output that MZmine rejects.
"""
import argparse
import hashlib
import re
from pathlib import Path


def rewrite_indexed_mzml(src, dst, fix=lambda line: line):
    """ Copy src -> dst passing every line through fix(), then rebuild the index and checksum """
    src, dst = Path(src), Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)

    offsets = {'spectrum': [], 'chromatogram': []}
    sha = hashlib.sha1()
    pos = 0
    with open(src, 'rb') as ifile, open(dst, 'wb') as ofile:
        def emit(text):
            nonlocal pos
            b = text.encode('utf-8')
            ofile.write(b)
            sha.update(b)
            pos += len(b)

        for raw in ifile:
            line = raw.decode('utf-8')
            stripped = line.lstrip()
            if stripped.startswith('<indexList'):
                break   # rebuilt below
            line = fix(line)
            stripped = line.lstrip()
            for kind in ('spectrum', 'chromatogram'):
                if stripped.startswith(f'<{kind} '):
                    ident = re.search(r' id="([^"]*)"', line).group(1)
                    offsets[kind].append((ident, pos + (len(line.encode('utf-8')) - len(stripped.encode('utf-8')))))
            emit(line)

        index_start = pos
        emit(f'  <indexList count="{sum(1 for v in offsets.values() if v)}">\n')
        for kind, entries in offsets.items():
            if not entries:
                continue
            emit(f'    <index name="{kind}">\n')
            for ident, off in entries:
                emit(f'      <offset idRef="{ident}">{off}</offset>\n')
            emit('    </index>\n')
        emit('  </indexList>\n')
        emit(f'  <indexListOffset>{index_start}</indexListOffset>\n')
        emit('  <fileChecksum>')
        ofile.write(f'{sha.hexdigest()}</fileChecksum>\n</indexedmzML>\n'.encode('utf-8'))


EMPTY_ARRAYS = """{i}<binaryDataArrayList count="2">
{i}  <binaryDataArray encodedLength="0">
{i}    <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value="" />
{i}    <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value="" />
{i}    <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value="" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z" />
{i}    <binary></binary>
{i}  </binaryDataArray>
{i}  <binaryDataArray encodedLength="0">
{i}    <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value="" />
{i}    <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value="" />
{i}    <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitCvRef="MS" unitAccession="MS:1000131" unitName="number of detector counts" />
{i}    <binary></binary>
{i}  </binaryDataArray>
{i}</binaryDataArrayList>
"""


def fixup_msconvert_mzml(src, dst):
    """
    Repairs msconvert output (seen converting Agilent mzXML):
    - <run> lacks the defaultInstrumentConfigurationRef that the mzML schema
      requires and MZmine's parser refuses to do without; point it at the
      file's first instrumentConfiguration
    - scans with no peaks are written with no binaryDataArrayList at all,
      which MassQL/pyteomics can't read; give them empty m/z and intensity
      arrays
    """
    instrument_id = None
    with open(src, encoding='utf-8') as ifile:
        for line in ifile:
            m = re.search(r'<instrumentConfiguration id="([^"]*)"', line)
            if m:
                instrument_id = m.group(1)
                break
            if '<run ' in line:
                break

    no_arrays = False   # inside a spectrum that hasn't shown a binaryDataArrayList yet

    def fix(line):
        nonlocal no_arrays
        stripped = line.lstrip()
        if instrument_id and stripped.startswith('<run ') and 'defaultInstrumentConfigurationRef' not in line:
            return line.replace('<run ', f'<run defaultInstrumentConfigurationRef="{instrument_id}" ', 1)
        if stripped.startswith('<spectrum '):
            no_arrays = True
        elif stripped.startswith('<binaryDataArrayList'):
            no_arrays = False
        elif stripped.startswith('</spectrum>') and no_arrays:
            no_arrays = False
            indent = line[:len(line) - len(stripped)]
            return EMPTY_ARRAYS.format(i=indent + '  ') + line
        return line

    rewrite_indexed_mzml(src, dst, fix)


def cli():
    argp = argparse.ArgumentParser(description=fixup_msconvert_mzml.__doc__)
    argp.add_argument('src')
    argp.add_argument('dst')
    args = argp.parse_args()
    fixup_msconvert_mzml(args.src, args.dst)


if __name__ == '__main__':
    cli()
