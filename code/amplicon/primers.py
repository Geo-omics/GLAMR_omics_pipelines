"""
Primer handling module
"""
import argparse
from dataclasses import asdict, dataclass, fields
import json
from pathlib import Path
from subprocess import Popen, PIPE
from typing import ClassVar


DEFAULT_AMPLICON_HMM_DB = '../data/reference/hmm_amplicons/combined.hmm'
DEFAULT_LOCATE_OUTPUT = 'primer_alignment_full.txt'
DEFAULT_JSON = 'primers.json'


def main():
    argp = argparse.ArgumentParser(description=__doc__)
    subp = argp.add_subparsers(dest='command', required=True)
    subp.add_parser(
        'print',
        help='Parse and print content of the google primer sheet',
    )
    updatep = subp.add_parser(
        'update',
        help='Write an updated primer_data.py file',
    )
    updatep.add_argument(
        '--out-json',
        default=DEFAULT_JSON,
        help='Json output file',
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
    locp.add_argument(
        '--output',
        default=DEFAULT_LOCATE_OUTPUT,
        help='Name of output file',
    )
    exgrp = argp.add_mutually_exclusive_group()
    exgrp.add_argument(
        '--sheet',
        help='The primer google sheet, saved as tsv',
    )
    exgrp.add_argument(
        '--json',
        help='Primers encoded as json, as made by the update command',
    )
    args = argp.parse_args()

    if args.sheet:
        primers = read_google_sheet(args.sheet)
    else:
        primers = Primer.load(args.json)

    match args.command:
        case 'print':
            print(*primers, sep='\n')
        case 'update':
            update_json(primers, output=args.out_json)
        case 'locate':
            locate_sequences(primers, db=args.hmm_db, output=args.output)


@dataclass
class Primer:
    name: str
    sequence: str
    direction: str
    gene_target: str
    citation: str
    hmm = None
    hmm_name: str
    start: int
    end: int
    region: str

    FWD: ClassVar = 'forward'
    REV: ClassVar = 'reverse'

    def __post_init__(self):
        if self.sequence:
            if not self.sequence.isupper():
                raise ValueError('sequence must be upper case only')
            if not self.sequence.isalpha():
                raise ValueError('sequence must be letters only')
            if not self.sequence.isascii():
                raise ValueError('sequence must be ascii only')

    def __str__(self):
        if self.direction == self.FWD:
            direc = '[fwd]'
        elif self.direction == self.REV:
            direc = '[rev]'
        else:
            direc = '<?>'

        if self.hmm_name:
            if self.start is None:
                start = '?'
            else:
                start = self.start
            if self.end is None:
                end = '?'
            else:
                end = self.end
            loc = f' {self.hmm_name}:{start}-{end}'
        else:
            loc = ''
        if loc and self.region:
            loc += f' ({self.region})'

        return f'{self.name} {direc}{loc}'

    @classmethod
    def load(cls, path=None):
        if path is None:
            path = Path(__file__).parent / DEFAULT_JSON
        with open(path) as ifile:
            return [Primer(**kw) for kw in json.load(ifile)]

    def is_forward(self):
        return self.direction == self.FWD

    def is_reverse(self):
        return self.direction == self.REV


def read_google_sheet(file_name):
    RENAMING = (
        ('primer_name', 'name'),
        ('primer_sequence', 'sequence'),
        ('HMM', 'hmm_name'),
    )
    field_names = [i.name for i in fields(Primer)]

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
            if not row:
                # google sheet does not save trailing empty cells, so empty
                # rows should have been skipped already, but let's be sure of
                # it
                continue
            row = row + [None] * (len(cols) - len(row))  # add missing fields
            row = {
                k: None if v == '' else v
                for k, v in zip(cols, row, strict=True)
            }

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

            if seq := row['sequence']:
                if seq in seqs:
                    print(f'WARNING: duplicate sequence: skipping line {lnum}')
                    continue
                else:
                    seqs.add(seq)
            else:
                row['sequence'] = None

            if row['start'] is not None:
                row['start'] = int(row['start'])

            if row['end'] is not None:
                row['end'] = int(row['end'])

            row = {k: v for k, v in row.items() if k in field_names}
            primers.append(Primer(**row))

    print(f'read {len(primers)} primer records')
    return primers


def update_json(primers, output=None):
    data = [asdict(i) for i in primers]
    data_txt = json.dumps(data, indent=4)
    if output:
        print(f'Writing to {output} ... ', end='', flush=True)
        with open(output, 'w') as ofile:
            ofile.write(data_txt)
        print('[Done]')
    else:
        print(data_txt)


def locate_sequences(primers, db=DEFAULT_AMPLICON_HMM_DB, output=None):
    if not Path(db).is_file():
        raise ValueError(f'hmm DB not found (no such file): {db}')
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
    if output is None:
        print(stdout.decode())
    else:
        with open(output, 'w') as ofile:
            ofile.write(stdout.decode())
        print(f'Output written to {output}')


if __name__ == '__main__':
    main()
