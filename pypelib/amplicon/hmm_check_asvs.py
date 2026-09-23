"""
Test ASV sequences against HMM DB
"""
import argparse
from collections import Counter
from itertools import groupby
from statistics import quantiles
from subprocess import run, DEVNULL
from tempfile import NamedTemporaryFile

from Bio import SeqIO

from .alignment import HMMRAlignment
from .hmm import HMM


def cli():
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument('hmm_db', help='combined hmm database')
    argp.add_argument('rep_seqs', help='Fasta file made with dada2 (the rep_seqs)')
    argp.add_argument('target-spec', help='The target spec')
    argp.add_argument(
        '--hmm-tblout',
        help='Optional name of (hmmr tabular formatted) alignment data output file',
    )
    argp.add_argument(
        '--final-asvs',
        help='Optional name of final ASV fasta file.',
    )
    args = argp.parse_args()
    main(args.hmm_db, args.fasta, args.target_spec, args.hmm_tblout, args.final_asvs)


def main(hmm_db_path, fasta_path, target_spec, aln_output, fasta_output,
         hmmr_threads=1):
    hmm, _, _ = HMM.spec2targets(target_spec)
    asvs = {i.id: i for i in SeqIO.parse(fasta_path, 'fasta')}

    print(f'Testing {len(asvs)} ASVs ...')

    with NamedTemporaryFile('rt') as tbl_out:
        cmd = [
            'nhmmscan',
            '--cpu', str(hmmr_threads - 1),  # worker threads
            '--tblout', tbl_out.name,
            hmm_db_path,
            fasta_path,
        ]
        print('running command:', *cmd)
        run(cmd, stdout=DEVNULL, check=True)
        alns0 = HMMRAlignment.load(tbl_out.name)

    alns0 = sorted(alns0, key=lambda x: x.qname)
    alns_out = []
    alns = {}
    for qname, grp in groupby(alns0, key=lambda x: x.qname):
        grp = list(grp)
        for aln in grp:
            if aln.model == hmm:
                grp.remove(aln)
                alns[qname] = aln
                alns_out.append(aln)
                break
        else:
            print(f'[WARNING] {qname}: no matching alignment! {grp}')
            alns_out += grp
            continue

        if others := [i for i in grp if i.model == hmm]:
            alns_out += others
            print(f'[WARNING] {qname}: multiple {hmm} alignments?! {others}')

    if all(i is None for i in alns.values()):
        raise RuntimeError(
            f'no rows with alignments to {hmm.name}, total lines: {len(alns0)}'
        )

    if nonaligned := [i for i in asvs if i not in alns]:
        print(f'[WARNING] {len(nonaligned)} ASVs w/o proper HMM alignment: ',
              ', '.join(nonaligned))

    if aln_output:
        HMMRAlignment.to_csv(alns_out, output=aln_output)
        print(f'[OK] {len(alns_out)} relevant alignments written to {aln_output}')

    attrs = ['hmmfrom', 'hmmto', 'envfrom', 'envto']
    data = {attname: [getattr(i, attname) for i in alns.values()] for attname in attrs}
    stats = zip(*[
        [
            min(nums),
            *[int(x) if x.is_integer() else x for x in quantiles(nums)],
            max(nums),
        ]
        for nums in data.values()
    ], strict=True)

    stat_header = ['[stats]', *attrs]
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

    strands = Counter(i.strand for i in alns.values())
    if len(strands) == 1:
        match list(strands.keys())[0]:
            case '+': pass
            case '-':
                print('[INFO] Sequences are on - strand, reverse-complementing ASVs!')
                for i in asvs:
                    asvs[i] = asvs[i].reverse_complement()
            case _: raise ValueError('invalid strand')
    else:
        raise RuntimeError(f'mixed strandedness: {strands}')

    with open(fasta_output, 'w') as ofile:
        for i in asvs:
            if aln := alns.get(i):
                envmin, envmax = sorted([aln.envfrom, aln.envto])
                seqid = (f'{i} {aln.model.name}:{aln.hmmfrom}..{aln.hmmto} '
                         f'env:{envmin}..{envmax}')
            else:
                seqid = f'{i} hmm_check_failed'

            ofile.write(f'>{seqid}\n{asvs[i].seq}\n')


if __name__ == '__main__':
    cli()
