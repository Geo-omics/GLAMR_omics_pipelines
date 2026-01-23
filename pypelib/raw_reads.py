"""
Check downloaded raw reads fastq files
"""
from contextlib import ExitStack
from functools import partial
from pathlib import Path
import subprocess

from .utils import load_stats


def check(params=None, stats=None, num_spots=None):
    """
    Run some checks on raw fastq files

    Returns None if everything went well.  Otherwise raises an exception.

    params:
        A params structure passed from snakemake.
    stats:
        tab-separated file screated with "seqkit stats -a -T".  If this is
        None, then is it assumed that params has it.
    """
    if stats is None:
        stats = params.stats_file
    if num_spots is None:
        num_spots = params.num_spots

    errs = []
    try:
        num_spots = int(num_spots)
    except (TypeError, ValueError) as e:
        raise RuntimeError(
            f'Failed accessing num_spots: {num_spots=} {e=}'
        ) from e

    for file, data in load_stats(stats).items():
        if (num_seqs := data['num_seqs']) != num_spots:
            errs.append(
                f'In file {file}: expected {num_spots} reads but got '
                f'{num_seqs}'
            )
    if errs:
        raise RuntimeError('\n'.join(errs))


def make_stats(input, output, keep_existing=False):
    """
    Generate some raw read statistics with "seqkit stats"

    input:
        Input filename or list of files.
    output:
        Output filename or snakemake file thingie.
    keep_existing [bool]:
        If True and the output file exists then do nothing.
    """
    output = Path(str(output))
    if keep_existing and output.is_file():
        return

    if isinstance(input, str):
        infiles = [input]
    else:
        # assume a list of snakemake input file thingies or similar
        infiles = [str(i) for i in input]

    cmd = ['seqkit', 'stats', '--quiet', '--basename', '-a', '-T', *infiles]

    with open(output, 'w+b') as ofile:
        subprocess.run(cmd, stdout=ofile, check=True)

        ofile.seek(0)

        for lnum, _ in enumerate(ofile, start=1):
            pass
        expect = len(infiles) + 1
        if lnum != expect:
            raise RuntimeError(
                f'[ERROR] seqkit stats wrote {lnum} lines but {expect} are '
                f'expected'
            )

    print(f'[OK] stats written to {output}')


def post_download(input, output, params, layout=None):
    if download_dir := input.get('download_dir'):
        download_dir = Path(download_dir)
    else:
        check(params)
        return

    out_files = [Path(i) for i in output]

    match layout:
        case 'PAIRED':
            fq_files = [
                download_dir / f'{params.srr_accession}_1.fastq',
                download_dir / f'{params.srr_accession}_2.fastq',
            ]
        case 'SINGLE':
            fq_files = [download_dir / f'{params.srr_accession}.fastq']
        case _:
            raise ValueError(f'invalid layout: {layout}')

    gz_files = [i.with_suffix('.fastq.gz') for i in fq_files]

    # compress if neeeded
    if all(i.is_file() for i in gz_files):
        pass
    elif all(i.is_file() for i in fq_files):
        subprocess.run(
            ['gzip', '--keep', '--fast', *(str(i) for i in fq_files)],
            check=True,
        )
    else:
        raise RuntimeError(
            f'[ERROR] Unexpected directory content, consider deleting and try '
            f'again: {download_dir}'
        )

    # Rename files:
    # This is so that the correct file names appears in the stats file, but
    # hard link as to keep the original files in place in case this operation
    # has to be re-done (mostly for testing?)
    out_files0 = [download_dir / i.name for i in out_files]
    for src, dest in zip(gz_files, out_files0, strict=True):
        dest.unlink(missing_ok=True)
        dest.hardlink_to(src)

    # the stats
    stats_file = Path(params.stats_file)
    stats_file0 = download_dir / stats_file.name
    make_stats(out_files0, stats_file0, keep_existing=False)
    check(params, stats=stats_file0)

    # move files into read directory
    for src, dest in zip(out_files0, out_files, strict=True):
        src.rename(dest)
    stats_file0.rename(stats_file)


post_download_paired = partial(post_download, layout='PAIRED')
post_download_single = partial(post_download, layout='SINGLE')
