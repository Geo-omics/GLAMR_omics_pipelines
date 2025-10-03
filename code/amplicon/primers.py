"""
Primer handling module
"""
import argparse
from dataclasses import dataclass
from pathlib import Path
from subprocess import Popen, PIPE

from .hmm import HMM


DEFAULT_AMPLICON_HMM_DB = '../data/reference/hmm_amplicons/combined.hmm'


def main():
    argp = argparse.ArgumentParser(description=__doc__)
    subp = argp.add_subparsers(dest='command', required=True)
    printp = subp.add_parser(
        'print',
        help='Parse and print content of the googlel. primer sheet',
    )
    locp = subp.add_parser(
        'locate',
        help='Run nhmmscan against combined HMM database to locate the primers'
    )
    locp.add_argument(
        '--hmm-db',
        default=DEFAULT_AMPLICON_HMM_DB,
        help='Path to the combined HMMR model database',
    )
    argp.add_argument(
        'google_primer_sheet',
        help='The primer google sheet, saved as tsv',
    )
    args = argp.parse_args()

    primers = read_google_sheet(args.google_primer_sheet)
    match args.command:
        case 'print':
            print(*primers, sep='\n')
        case 'locate':
            locate_sequences(primers, db=args.hmm_db)


@dataclass
class Primer:
    name: str
    sequence: str
    direction: str
    gene_target: str
    citation: str
    hmm: HMM
    start: int
    end: int
    region: str

    def __post_init__(self):
        if self.sequence:
            if not self.sequence.isupper():
                raise ValueError('sequence must be upper case only')
            if not self.sequence.isalpha():
                raise ValueError('sequence must be letters only')
            if not self.sequence.isascii():
                raise ValueError('sequence must be ascii only')

    def __str__(self):
        if self.direction == 'forward':
            direc = '[fwd]'
        elif self.direction == 'reverse':
            direc = '[rev]'
        else:
            direc = '<?>'

        if self.hmm:
            if self.start is None:
                start = '?'
            else:
                start = self.start
            if self.end is None:
                end = '?'
            else:
                end = self.end
            loc = f' {self.hmm}:{start}-{end}'
        else:
            loc = ''
        if loc and self.region:
            loc += f' ({self.region})'

        return f'{self.name} {direc}{loc}'


def read_google_sheet(file_name):
    RENAMING = (
        ('primer_name', 'name'),
        ('primer_sequence', 'sequence'),
        ('HMM', 'hmm'),
    )

    primers = []
    with open(file_name) as ifile:
        cols = ifile.readline().strip().split('\t')

        names = set()
        seqs = set()
        for lnum, line in enumerate(ifile, start=2):
            line = line.strip()
            if not line or line.startswith('#'):
                # ignore empty or comment lines
                continue
            row = line.split('\t')
            row = row + [None] * (len(cols) - len(row))  # add missing fields
            row = {k: v for k, v in zip(cols, row, strict=True)}

            for old, new in RENAMING:
                row[new] = row.pop(old)

            name = row['name']
            if not name:
                print(f'WARNING: name missing on line {lnum}, skipping')
                continue
            if name in names:
                print(f'WARNING: duplicate name: {name} (skipping line '
                      f'{lnum})')
                continue
            else:
                names.add(name)

            seq = row['sequence']
            if seq:
                if seq in seqs:
                    print(f'WARNING: duplicate sequence: skipping line {lnum}')
                    continue
                else:
                    seqs.add(seq)

            if primer_pos := row.pop('primer_pos'):
                start, _, end = primer_pos.partition('-')
                try:
                    row['start'] = int(start)
                    row['end'] = int(end)
                except Exception as e:
                    print(f'WARNING: line {lnum}: failed parsing primer '
                          f'positions: {e}')
                    row['start'] = row['end'] = None

            primers.append(Primer(**row))

    print(f'read {len(primers)} primer records')
    return primers


def locate_sequences(primers, db=DEFAULT_AMPLICON_HMM_DB):
    if Path(db).is_file():
        raise ValueError('no such file:')
    data = b''.join(f'>{i.name}\n{i.sequence}\n'.encode() for i in primers)
    cmd = [
        'nhmmscan', '--cpu', '4', '--tblout', 'primer_alignments.csv',
        '--max',  # !important
        db, '-',
    ]
    p = Popen(cmd, stdin=PIPE, stdout=PIPE)
    stdout, stderr = p.communicate(data)
    if stderr:
        print(f'ERROR: {stderr.decode()}')
    with open('primer_alignment_full.txt', 'w') as ofile:
        ofile.write(stdout.decode())


if __name__ == '__main__':
    main()
