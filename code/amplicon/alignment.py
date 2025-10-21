from dataclasses import dataclass
from typing import ClassVar

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

    @property
    def pass_test(self):
        """ Tell if this this a good alignment """
        return self.score >= self.MIN_SCORE

    def primer_score(self):
        best_outer_names = None
        best_outer_score = None
        best_inner_names = None
        best_inner_score = None

        if self.strand == '+' and self.envfrom < self.envto:
            # alignment based on forward read
            forward = True
            primers = self.model.fwd_primers
        elif self.strand == '-' and self.envto < self.envfrom:
            # alignment based on reverse read
            forward = False
            primers = self.model.rev_primers
        else:
            raise RuntimeError('unclear directionality')

        self.direction = 'fwd' if forward else 'rev'

        if not primers:
            self.has_primer = None
            self.nearest_primer = None
            self.primer_score = None
            return

        for i in primers:
            # scores are positive if read goes beyond the primer
            # score is negative if read is short of primer
            # score of zero is a perfect alignment
            if forward:
                outer_score = i.start - self.hmmfrom
                inner_score = i.end - self.hmmfrom
            else:
                outer_score = self.hmmto - i.end
                inner_score = self.hmmto - i.start

            if best_outer_names is None:
                best_outer_names = [i.name]
                best_inner_names = [i.name]
                best_outer_score = outer_score
                best_inner_score = inner_score
            else:
                if abs(outer_score) < abs(best_outer_score):
                    best_outer_score = outer_score
                    best_outer_names = [i.name]
                elif abs(outer_score) == abs(best_outer_score):
                    best_outer_names.append(i.name)

                if abs(inner_score) < abs(best_inner_score):
                    best_inner_score = inner_score
                    best_inner_names = [i.name]
                elif abs(inner_score) == abs(best_inner_score):
                    best_inner_names.append(i.name)

        if abs(best_outer_score) <= abs(best_inner_score):
            has_primer = True
            names = ','.join(sorted(best_outer_names))
            score = best_outer_score
        else:
            has_primer = False
            names = ','.join(sorted(best_inner_names))
            score = best_inner_score

        self.has_primer = has_primer
        self.nearest_primer = names
        self.primer_score = score
