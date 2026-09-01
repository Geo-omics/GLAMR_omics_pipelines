"""
Detect the vendor format of a raw LC-MS/MS file (or directory) and decide
how it needs to be converted to mzML.

New vendor formats are added here, not scattered across the Snakefile.
"""
import argparse
from dataclasses import dataclass
from pathlib import Path

from ..utils import UsageError


THERMO_MAGIC = 'Finnigan'.encode('utf-16-le')

# Directory-format vendor markers: (path glob relative to the .d dir) -> vendor
BRUKER_D_MARKERS = ('analysis.tdf', 'analysis.tsf', 'analysis.baf')
AGILENT_D_MARKERS = ('AcqData', 'MSScan.bin', 'MSPeak.bin')

# Formats that are known to exist but that our conversion route
# (ProteoWizard under Wine, docker://chambm/pwiz-skyline-i-agree-to-the-vendor-licenses)
# does not support. Fail loudly rather than silently emitting an empty mzML.
KNOWN_UNSUPPORTED = {
    '.t2d': 'ABI/SCIEX T2D is not supported by our msconvert route',
    '.yep': 'Bruker YEP (bare, not in a .d directory) is not supported by our msconvert route',
    '.fid': 'Bruker FID (bare, not in a .d directory) is not supported by our msconvert route',
}

# formats mzmine/metaMS can read directly, no conversion needed
PASSTHROUGH_SUFFIXES = {'.mzml', '.mzxml'}


@dataclass(frozen=True)
class FormatInfo:
    vendor: str
    """ short vendor/format tag, e.g. 'thermo_raw', 'bruker_d', 'mzml' """
    source_path: Path
    needs_conversion: bool
    conversion_tool: str | None
    """ None (no conversion needed), 'thermorawfileparser', or 'msconvert' """


def cli():
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument('path', help='Path to a raw spectra file or directory')
    args = argp.parse_args()
    try:
        info = detect(args.path)
    except UsageError as e:
        argp.error(e)
    else:
        print(f'vendor: {info.vendor}')
        print(f'needs_conversion: {info.needs_conversion}')
        print(f'conversion_tool: {info.conversion_tool}')


def detect(path) -> FormatInfo:
    """ Identify the vendor format of a raw spectra file or directory """
    path = Path(path)
    if not path.exists():
        raise UsageError(f'no such file or directory: {path}')

    if path.is_dir():
        return _detect_directory_format(path)
    return _detect_file_format(path)


def _detect_directory_format(path: Path) -> FormatInfo:
    suffix = path.suffix.lower()
    if suffix == '.d':
        if any((path / m).exists() for m in BRUKER_D_MARKERS):
            # mzmine reads Bruker .d natively, but we still convert to mzML
            # for uniformity with metaMS and with harvest_metadata's mzML
            # header parsing.
            return FormatInfo('bruker_d', path, True, 'msconvert')
        if any((path / m).exists() for m in AGILENT_D_MARKERS):
            return FormatInfo('agilent_d', path, True, 'msconvert')
        raise UsageError(
            f'unrecognized .d directory format (neither Bruker nor Agilent '
            f'markers found): {path}'
        )
    if suffix == '.raw':
        # Waters raw data is a directory, unlike Thermo's single-file .raw
        return FormatInfo('waters_raw', path, True, 'msconvert')
    raise UsageError(f'unrecognized directory format: {path}')


def _detect_file_format(path: Path) -> FormatInfo:
    suffix = path.suffix.lower()

    if suffix in KNOWN_UNSUPPORTED:
        raise UsageError(
            f'{path}: {KNOWN_UNSUPPORTED[suffix]}. This file cannot be '
            f'processed by this pipeline.'
        )

    if suffix in PASSTHROUGH_SUFFIXES:
        return FormatInfo(suffix.removeprefix('.'), path, False, None)

    if suffix == '.raw':
        with open(path, 'rb') as ifile:
            head = ifile.read(len(THERMO_MAGIC) + 16)
        if THERMO_MAGIC not in head:
            raise UsageError(
                f'{path}: has a .raw extension but does not look like a '
                f'Thermo raw file (missing Finnigan magic bytes)'
            )
        return FormatInfo('thermo_raw', path, True, 'thermorawfileparser')

    if suffix in ('.wiff', '.wiff2'):
        return FormatInfo('sciex_wiff', path, True, 'msconvert')

    if suffix == '.lcd':
        return FormatInfo('shimadzu_lcd', path, True, 'msconvert')

    raise UsageError(f'unrecognized file format: {path} (suffix {suffix!r})')


if __name__ == '__main__':
    cli()
