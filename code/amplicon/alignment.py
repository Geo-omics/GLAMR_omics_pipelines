from dataclasses import dataclass
from functools import cached_property
from typing import ClassVar

from . import get_models
from .hmm import HMM


class MissingPrimerInfo(Exception):
    pass


@dataclass
class HMMRAlignment:
    """ An nhmmscan tab output row """
    model: HMM
    qname: str
    hmmfrom: int
    hmmto: int
    envfrom: int
    envto: int
    strand: str
    score: float

    MIN_SCORE: ClassVar[float] = 10.0

    def __post_init__(self):
        if self.strand == '+' and self.envfrom < self.envto:
            # alignment based on forward read
            self.direction = 'fwd'
        elif self.strand == '-' and self.envto < self.envfrom:
            # alignment based on reverse read
            self.direction = 'rev'
        else:
            raise RuntimeError('unclear directionality')

    @property
    def pass_test(self):
        """ Tell if this this a good alignment """
        return self.score >= self.MIN_SCORE

    def score_fwd_primers(self):
        return sorted(PrimerMatch(self, i) for i in self.model.fwd_primers)

    def score_rev_primers(self):
        return sorted(PrimerMatch(self, i) for i in self.model.rev_primers)

    @cached_property
    def fwd_match(self):
        """ Get the top-scoring match """
        matches = self.score_fwd_primers()
        if matches:
            return matches[0]
        else:
            return None

    @cached_property
    def rev_match(self):
        """ Get highest-scoring rev primer match to the right of fwd match """
        # TODO: if the best fwd_match is much worse than the best rev_match we
        # but the might want to consider throwing out the fwd_match
        matches = self.score_rev_primers()
        if self.fwd_match is None:
            if matches:
                return matches[0]
            else:
                return None

        for i in matches:
            # only consider matches to right of fwd primer
            if self.fwd_match.primer.end < i.primer.start:
                return i
        return None

    @classmethod
    def load(cls, path):
        """
        Load alignmet data from a HMMR tblout formatted file

        :rtype: list of HMMRAlignment instances
        """
        models = get_models()
        rows = []
        with open(path) as ifile:
            for lnum, line in enumerate(ifile, start=1):
                if line.startswith('#'):
                    # two header lines and some trailing comments
                    continue
                row = line.strip().split()
                try:
                    model_name = row[0]
                except IndexError as e:
                    raise RuntimeError(f'empty line at {lnum}') from e

                try:
                    model = models[model_name]
                except KeyError as e:
                    raise RuntimeError(
                        f'unknown HMM model: "{model_name}" at line {lnum}'
                    ) from e

                try:
                    row = cls(
                        model=model,
                        qname=row[2],
                        hmmfrom=int(row[4]),
                        hmmto=int(row[5]),
                        envfrom=int(row[8]),
                        envto=int(row[9]),
                        strand=row[11],
                        score=float(row[13]),
                    )
                except (IndexError, ValueError) as e:
                    raise RuntimeError(
                        f'failed parsing line {lnum}: {e} -- \n{row=}'
                    ) from e
                rows.append(row)
        return rows


class PrimerMatch:
    def __init__(self, alignment, primer):
        if primer.hmm is not alignment.model:
            raise ValueError('alignment vs. primer model mismatch')

        self.primer = primer

        # is primer ahead or trailing?
        is_ahead = primer.end < alignment.hmmfrom
        is_trailing = alignment.hmmto < primer.start
        if is_ahead and is_trailing:
            raise RuntimeError('cannot be both: bad coordinates?')
        self.clean = is_ahead or is_trailing

        # The score: generally, it's the offset of the read w.r.t the primer.
        # A zero score means the read lines up exactly with the primer (but in
        # two different senses.)
        if primer.is_forward():
            if is_ahead:
                # score is <=0
                self.score = alignment.hmmfrom - primer.end - 1
            elif is_trailing:
                # score is >= alignment length
                self.score = alignment.hmmfrom - primer.start
            else:
                # dirty / score is offset
                self.score = alignment.hmmfrom - primer.start
        elif primer.is_reverse():
            if is_ahead:
                # score is <0 smaller than negative alignment length
                self.score = alignment.hmmto - primer.end
            elif is_trailing:
                # score is >=0
                self.score = alignment.hmmto - primer.start + 1
            else:
                # dirty
                self.score = alignment.hmmto - primer.end
        else:
            raise RuntimeError('logic error')

    def __str__(self):
        state = 'clean' if self.clean else 'dirty'
        sign = '+' if 0 <= self.score else ''
        return f'{self.primer.name}/{state}{sign}{self.score}'

    def __eq__(self, other):
        return abs(self.score) == abs(other.score)

    def __lt__(self, other):
        # less is "better" and zero is best score
        return abs(self.score) < abs(other.score)
