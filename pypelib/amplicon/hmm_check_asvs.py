"""
Test ASV sequences against HMM DB
"""
import argparse
from statistics import quantiles
from subprocess import run, DEVNULL
from tempfile import NamedTemporaryFile
import sys

from .alignment import HMMRAlignment
from .hmm import HMM


def cli():
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument('hmm_db', help='combined hmm database')
    argp.add_argument('target_spec', help='e.g. 16S_rRNA_bac.524-776')
    argp.add_argument('asv_fasta', help='Fasta file with ASV sequences')
    argp.add_argument('--output', help='Output file')
    args = argp.parse_args()
    main(args.hmm_db, args.target_spec, args.asv_fasta, output=args.output)


def main(hmm_db_path, target_spec, fasta_path, output=None):
    with open(fasta_path) as ifile:
        fasta_count = sum(1 for line in ifile if line.startswith('>'))
    print(f'Testing {fasta_count} ASVs ...')

    hmm, _, _ = HMM.spec2targets(target_spec)
    with NamedTemporaryFile('rt') as tbl_out:
        cmd = [
            'nhmmscan',
            '--tblout', tbl_out.name,
            '--cpu', '2',
            hmm_db_path,
            fasta_path,
        ]
        run(cmd, stdout=DEVNULL, check=True)
        all_alns = HMMRAlignment.load(tbl_out.name)
        alns = [i for i in all_alns if i.model == hmm]
        if alns:
            print(f'Found {len(alns)} matching ASV alignments')
        else:
            raise RuntimeError(
                f'no matching alignments at all out of a total of {len(all_alns)}'
            )

        cols = ['qname', 'hmmfrom', 'hmmto', 'envfrom', 'envto', 'strand', 'score']
        ofile = open(output, 'w') if output else sys.stdout
        ofile.write('\t'.join(cols))
        ofile.write('\n')
        for i in alns:
            row = (str(getattr(i, attr)) for attr in cols)
            ofile.write('\t'.join(row))
            ofile.write('\n')
        if output:
            ofile.close()
            print(f'[OK] {output} written')

        """

        print(f'BOEK {len(alns)=}')
        header = tbl_out.readline().rstrip('\n')
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
        """

    stats_attrs = ('hmmfrom', 'hmmto', 'envfrom', 'envto')
    stats = zip(*[
        [
            min(nums),
            *[int(x) if x.is_integer() else x for x in quantiles(nums)],
            max(nums),
        ]
        for nums in ([getattr(al, attr) for al in alns] for attr in stats_attrs)
    ], strict=True)

    stat_header = ['[stats]', *stats_attrs]
    names = ['min', 'Q1', 'Q2', 'Q3', 'max']
    widths = [len(i) for i in stat_header]
    print(*stat_header, sep=' ')
    for name, row in zip(names, stats):
        row = [
            str(val).rjust(width)
            for val, width in zip([name, *row], widths)
        ]
        print(*row, sep=' ')
    print('total:', len(alns))


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
