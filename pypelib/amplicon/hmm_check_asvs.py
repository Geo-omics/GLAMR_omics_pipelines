"""
Test ASV sequences against HMM DB
"""
import argparse
from statistics import quantiles
from subprocess import run, DEVNULL
from tempfile import NamedTemporaryFile

from .hmm import HMM


def cli():
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument('hmm_db', help='combined hmm database')
    argp.add_argument('asv_fasta', help='Fasta file with ASV sequences')
    args = argp.parse_args()
    main(args.hmm_db, args.fasta)


def main(hmm_db_path, fasta_path, target_spec, output=None):
    with open(fasta_path) as ifile:
        fasta_count = sum(1 for line in ifile if line.startswith('>'))
    print(f'Testing {fasta_count} ASVs ...')

    colnames = ['hmmfrom', 'hmm to', 'envfrom', 'env to']
    hmm, _, _ = HMM.spec2targets(target_spec)
    data = {col: [] for col in colnames}
    with NamedTemporaryFile('rt') as tbl_out, open(output, 'w') as ofile:
        cmd = [
            'nhmmscan',
            '--tblout', tbl_out.name,
            hmm_db_path,
            fasta_path,
        ]
        run(cmd, stdout=DEVNULL, check=True)
        header = tbl_out.readline().rstrip('\n')
        ofile.write(header + '\n')
        slices = get_cols(colnames, header)
        for lnum, line in enumerate(tbl_out, start=2):
            line = line.rstrip()
            if not line.startswith(hmm.name):
                continue

            for col in colnames:
                value = line[slices[col]]
                try:
                    numval = int(value.strip())
                except ValueError as e:
                    print(f'[Error] {col=} {value=}\n{header}\n{line}')
                    raise ValueError(f'on line {lnum}: {e}') from e
                data[col].append(numval)
                ofile.write(line + '\n')

    print(f'Found {len(data[colnames[0]])} matching ASV alignments')
    print(f'[OK] {output} written')

    if not data[colnames[0]]:
        raise RuntimeError(
            f'no rows with alignments to {hmm.name}, total lines: {lnum}'
        )

    stats = zip(*[
        [
            min(nums),
            *[int(x) if x.is_integer() else x for x in quantiles(nums)],
            max(nums),
        ]
        for nums in data.values()
    ], strict=True)

    stat_header = ['[stats]', *colnames]
    names = ['min', 'Q1', 'Q2', 'Q3', 'max']
    widths = [len(i) for i in stat_header]
    print(*stat_header, sep=' ')
    for name, row in zip(names, stats):
        row = [
            str(val).rjust(width)
            for val, width in zip([name, *row], widths)
        ]
        print(*row, sep=' ')
    print('total:', len(data[colnames[0]]))


def get_cols(colnames, headerline):
    """
    Figure out column coordinates for HMMR tabular output file format
    """
    cols = {}
    for col in colnames:
        try:
            start = headerline.index(col)
        except ValueError as e:
            raise ValueError(f'not a column name: {col} ({e})') from e
        # find stop: the beginning of next column
        stop = start + len(col)
        while stop < len(headerline) and headerline[stop] == ' ':
            stop += 1
        cols[col] = slice(start, stop)
    return cols


if __name__ == '__main__':
    cli()
