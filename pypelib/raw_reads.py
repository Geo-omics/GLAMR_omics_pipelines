"""
Check downloaded raw reads fastq files
"""
from contextlib import ExitStack
from functools import partial
import gzip
from itertools import batched
import json
from pathlib import Path
import re
import shutil
import subprocess

from .utils import load_stats, logme


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


def parse_runinfo(file, key=None):
    """
    Get value from runinfo file

    Helper for rules that need bits from the runinfo file.
    The runinfo file is a json-formatted file obtained by "kingfisher annotate"

    If key is None, then all data is returned.
    """
    with ExitStack() as stack:
        try:
            ifile = open(file)
        except TypeError:
            # file is file-like, e.g. TextIOWrapper, got here via parse_input()
            ifile = file
            ifile.seek(0)  # file handle re-used across multiple calls
        else:
            # file is Path or snakemake _IOFile, e.g. got here via input
            # function
            stack.enter_context(ifile)

        data = json.load(ifile)  # expecting a list
        if len(data) == 0:
            raise RuntimeError(f'empty runinfo? {ifile.name}')
        if len(data) > 1:
            # To be implemented if seen in the wild
            raise RuntimeError(f'multiple exp packages? runinfo: {ifile.name}')
        data = data[0]  # should get a dict

        if key is None:
            return data
        else:
            return data[key]


@logme()
def post_download(input, output, params, layout=None):
    if download_dir := input.get('download_dir'):
        download_dir = Path(download_dir)
    else:
        check(params)
        return

    out_files = [Path(i) for i in output]

    fwd_fq = download_dir / f'{params.srr_accession}_1.fastq'
    rev_fq = download_dir / f'{params.srr_accession}_2.fastq'
    sgl_fq = download_dir / f'{params.srr_accession}.fastq'

    match layout:
        case 'PAIRED':
            fq_files = [fwd_fq, rev_fq]
        case 'SINGLE':
            fq_files = [sgl_fq]
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
            f'[ERROR] ({params=} {layout=} Unexpected directory content, '
            f'consider deleting and try again: {download_dir}'
        )

    num_spots = params.num_spots
    print_interleave_msg = False

    if layout == 'SINGLE':
        try:
            inter_fq_gz = gz_files[0]
            check_interleaved(inter_fq_gz)
        except NotInterleaved:
            # proceed normally
            pass
        else:
            # Looks like SRA meta data got it wrong - it's paired
            gz_files = [
                fwd_fq.with_suffix('.fastq.gz'),
                rev_fq.with_suffix('.fastq.gz'),
            ]
            print('De-interleaving... ', end='', flush=True)
            deinterleave(inter_fq_gz, *gz_files)
            print('[OK]')
            # save as paired output files, even though those were not requested
            orig_out = out_files[0]
            if 'single' not in orig_out.name:
                raise RuntimeError(
                    f'do not know how to rename single raw fastq file into fwd'
                    f' and rev file pair: {out_files=}'
                )
            fwd_out_name = orig_out.name.replace('single', 'fwd')
            rev_out_name = orig_out.name.replace('single', 'rev')
            out_files = [orig_out.with_name(fwd_out_name),
                         orig_out.with_name(rev_out_name)]

            # write a new minimal runinfo.json, keep a backup
            runinfo = Path(input.get('runinfo'))
            orig_info = parse_runinfo(runinfo)
            num_spots = orig_info['spots'] // 2
            info = {
                'library_layout': 'PAIRED',
                'run': orig_info['run'],
                'spots': num_spots,
            }
            backup = runinfo.with_suffix(runinfo.suffix + '.orig')
            shutil.copy2(runinfo, backup)
            with open(runinfo, 'w') as ofile:
                json.dump([info], ofile, indent=4)
                ofile.write('\n')
            print(f'Fixed {runinfo} -- kept backup under {backup}')
            print_interleave_msg = True

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
    check(params, stats=stats_file0, num_spots=num_spots)

    # move files into read directory
    for src, dest in zip(out_files0, out_files, strict=True):
        src.rename(dest)
    stats_file0.rename(stats_file)

    if print_interleave_msg:
        print(
            f'[NOTICE] All went well, but rule for output {output} is going to'
            f' fail since a de-interleaved pair of raw fastq files got made '
            f'and not the expected single raw fastq.  Simply re-run your '
            f'downstream rule to continue.  You may safely delete '
            f'{download_dir}'
        )


post_download_paired = partial(post_download, layout='PAIRED')
post_download_single = partial(post_download, layout='SINGLE')


class NotInterleaved(Exception):
    pass


fastq_id_pat = re.compile(
    r'^@(?P<srr_accn>\w+)\.(?P<seqnum>[0-9]+) (?P<pair_id>.*)$'
)
"""
Raw fastq sequence ID pattern, e.g.:
@SRR10522209.1 MISEQ04_297_000000000-AMFNM_1_1101_18742_1259 length=301
or
@SRR10522209.1 MISEQ04_297_000000000-AMFNM_1_1101_15772_1034/1
"""


def check_interleaved(fastq_file):
    """
    Check if given fastq file is in interleaved format

    Sequence indentifier lines must match exactly and be unique for each pair
    except for the sequential number.

    Returns number of read pairs if the file is in interleaved format, raises
    NotInterleaved otherwise.
    """
    fastq_file = Path(fastq_file)
    with ExitStack() as estack:
        if fastq_file.suffix == '.gz':
            ifile = gzip.open(fastq_file, mode='rt')
        else:
            ifile = open(fastq_file)
        estack.enter_context(ifile)

        pair_ids = set()
        accn = None
        records = batched(ifile, n=8, strict=True)
        for recn, lines in enumerate(records, start=1):
            id1, seq1, plus1, qual1, id2, seq2, plus2, qual2 = lines

            m1 = m2 = 'n/a'

            def errinfo():
                return (f'at line {recn * 8} in {fastq_file}\n'
                        f'Match 1 {m1=}\nMatch 2 {m2=}')

            m1 = fastq_id_pat.match(id1)
            if m1 is None:
                raise NotInterleaved(
                    f'failed parsing sequence identifier: {id1=} - {errinfo()}'
                )
            m1 = m1.groupdict()

            m2 = fastq_id_pat.match(id2)
            if m2 is None:
                raise NotInterleaved(
                    f'failed parsing sequence identifier: {id2=} - {errinfo()}'
                )
            m2 = m2.groupdict()

            seqnum1 = int(m1['seqnum'])
            seqnum2 = int(m2['seqnum'])
            if not seqnum1 + 1 == seqnum2:
                raise NotInterleaved(
                    f'non-consequitive sequence numbers - {errinfo()}'
                )

            if accn is None:
                accn = m1['srr_accn']
            if not (accn == m1['srr_accn'] == m2['srr_accn']):
                raise NotInterleaved(f'accession mismatch -- {errinfo()}')

            if m1['pair_id'] != m2['pair_id']:
                raise NotInterleaved(f'pair IDs mismatch -- {errinfo()}')

            if m1['pair_id'] in pair_ids:
                raise NotInterleaved(f'duplicate pair ID -- {errinfo()}')

            pair_ids.add(m1['pair_id'])

    if not pair_ids:
        raise RuntimeError('empty input file?')
    print(f'Fastq file {fastq_file} is interleaved with {recn} read pairs')
    return recn


def fix_reverse_seq_id(fastq_seq_id):
    """
    subtract one from sequence number

    Helper function for deinterleave() to fix fastq headers.
    """
    m = fastq_id_pat.match(fastq_seq_id)
    if m is None:
        raise ValueError(f'error parsing fastq seq ID line: {fastq_seq_id}')

    srr_accn, seqnum, pair_id = m.groups()
    seqnum = int(seqnum) - 1  # previous one
    return f'@{srr_accn}.{seqnum} {pair_id}\n'


def deinterleave(fastq_in, fastq_out_fwd, fastq_out_rev):
    """
    De-interleaves a fastq file

    The process is rather blind, goes by sets of four lines, just copies them
    to the output files.  Input file must be gzipped.
    """
    with gzip.open(fastq_in, 'rt') as ifile:
        with gzip.open(fastq_out_fwd, 'wt') as ofwd:
            with gzip.open(fastq_out_rev, 'wt') as orev:
                for lines in batched(ifile, n=8):
                    lines = list(lines)
                    for i in lines[:4]:
                        ofwd.write(i)
                    lines[4] = fix_reverse_seq_id(lines[4])
                    for i in lines[4:]:
                        orev.write(i)
