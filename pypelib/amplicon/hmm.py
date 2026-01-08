import re


class HMM:
    """
    HMM model
    """
    def __init__(self, name, regions=None, primers=None, length=None):
        if not name:
            raise ValueError('name must be supplied')

        if not isinstance(name, str):
            raise TypeError('name must be a str')

        self.name = name
        if regions is None:
            regions = {}

        # check regions consistency
        last = 0
        for start, end in regions.values():
            if start <= 0 or end <= 0:
                raise ValueError('coords must be positive')
            if start <= last:
                raise ValueError(
                    'start must be larger than previous end (and >0)'
                )
            if end < start:
                raise ValueError('end must be larger than start')
            if length is not None and length < end:
                raise ValueError('end must not be larger than length')
            last = end
        self.regions = regions
        self.length = length

        fwd_primers = []
        rev_primers = []
        if primers is not None:
            for i in primers:
                if i.start is None or i.end is None:
                    continue
                if length is not None and length < i.end:
                    raise ValueError('prim end must not be larger than length')
                if i.hmm is not None and i.hmm.name != name:
                    raise ValueError(
                        'a different HMM is already assigned to this primer'
                    )
                else:
                    i.hmm = self
                if i.is_forward():
                    fwd_primers.append(i)
                elif i.is_reverse():
                    rev_primers.append(i)
                else:
                    # skip
                    continue
        self._fwd_primers = tuple(fwd_primers)
        self._rev_primers = tuple(rev_primers)
        self._fwd_primer = {i.name: i for i in fwd_primers}
        self._rev_primer = {i.name: i for i in rev_primers}

    def __str__(self):
        return self.name

    def __repr__(self):
        return f'<{self.__class__.__name__}: {self.name}>'

    @property
    def fwd_primers(self):
        return self._fwd_primers

    @property
    def rev_primers(self):
        return self._rev_primers

    @property
    def fwd_primer(self):
        """ provide access to primers by name """
        return self._fwd_primer

    @property
    def rev_primer(self):
        """ provide access to primers by name """
        return self._rev_primer

    def format_target_spec(self, fwd_primer, rev_primer):
        spec = self.name
        if fwd_primer is not None and rev_primer is not None:
            start, end = self.get_target_interval(fwd_primer, rev_primer)
            spec += f'.{start}-{end}'
        return spec

    def format_target(self, fwd_primer, rev_primer):
        if isinstance(fwd_primer, str):
            fwd_primer = self.fwd_primer[fwd_primer]
        elif fwd_primer not in self.fwd_primers:
            raise ValueError(f'not a fwd primer for {self}')
        if isinstance(rev_primer, str):
            rev_primer = self.rev_primer[rev_primer]
        elif rev_primer not in self.rev_primers:
            raise ValueError(f'not a rev primer for {self}')

        return f'{self.name}.{fwd_primer.name}.{rev_primer.name}'

    @classmethod
    def parse_target(cls, target):
        """
        Return hmm and primer instances corresponding to target string

        This is the inverse of format_target().
        """
        from . import get_models
        hmm_name, *primers = target.split('.')
        fwd_primer, rev_primer = primers if primers else (None, None)
        try:
            obj = get_models()[hmm_name]
        except KeyError as e:
            raise LookupError(f'Invalid HMM name: {e}') from e

        if fwd_primer:
            try:
                fwd_primer = obj.fwd_primer[fwd_primer]
            except KeyError as e:
                raise LookupError(f'not a fwd primer for {obj}: {e}') from e

        if rev_primer:
            try:
                rev_primer = obj.rev_primer[rev_primer]
            except KeyError as e:
                raise LookupError(f'not a rev primer for {obj}: {e}') from e
        return obj, fwd_primer, rev_primer

    def get_target_interval(self, fwd_primer, rev_primer):
        """ Get start and end positions of the exclusive sequence targeted by
        given primer pair """
        if isinstance(fwd_primer, str):
            fwd_primer = self.fwd_primer[fwd_primer]
        if isinstance(rev_primer, str):
            rev_primer = self.rev_primer[rev_primer]
        return fwd_primer.end + 1, rev_primer.start - 1

    @classmethod
    def spec2targets(cls, spec):
        """ Get list of matching targets from given target spec """
        from . import get_models

        pat = re.compile(
            r'^(?P<name>\w+)(.(?P<start>[0-9]+)-(?P<end>[0-9]+))?$'
        )
        m = pat.match(spec)
        if m is None:
            raise ValueError(
                f'Does not conform to target spec syntax: "{spec}"'
            )
        m = m.groupdict()

        hmm_name = m['name']
        start = m['start']
        end = m['end']

        try:
            hmm = get_models()[hmm_name]
        except KeyError as e:
            raise ValueError(f'Invalid HMM name: {e}') from e

        if start is None and end is None:
            return [hmm.name]

        start = int(start)
        end = int(end)

        fwd_p_names = [i.name for i in hmm.fwd_primers if i.end + 1 == start]
        rev_p_names = [i.name for i in hmm.rev_primers if i.start - 1 == end]

        if fwd_p_names and rev_p_names:
            return [
                f'{hmm_name}.{i}.{j}'
                for i in sorted(fwd_p_names)
                for j in sorted(rev_p_names)
            ]
        else:
            msg = []
            if not fwd_p_names:
                msg.append('No matching forward primers found.')
            if not rev_p_names:
                msg.append('No matching reverse primers found.')
            print(f'{hmm=}')
            for i in hmm.fwd_primers:
                print(f'  fwd: {i} end={i.end}')
            for i in hmm.rev_primers:
                print(f'  rev: {i} start={i.start}')
            raise ValueError(' '.join(msg))

    def get_regions(self, qstart, qend):
        """ get regions contained in interval -- any overlap counts """
        if qstart <= 0 or qend <= 0:
            raise ValueError('must be positive')
        if qend <= qstart:
            raise ValueError('qstart must be smaller than qend')
        if self.length is not None and self.length < qend:
            raise ValueError('qend must be smaller then model length')

        regions = []
        found_all = False
        for name, (start, end) in self.regions.items():
            if qend < start or end < qstart:
                # total miss
                if regions:
                    found_all = True
            else:
                # some overlap
                if found_all:
                    # assuming regions are monotonic, sanity check
                    raise RuntimeError('logic error')
                regions.append(name)
        return tuple(regions)

    def contains_primers(self, qstart, qend):
        """
        tell primer content of interval

        Interval start or end must be within 5 positions of primer.
        """
        fwd = [
            i for i in self.fwd_primers
            if (i.start - 5) <= qstart <= (i.start + 5) and qstart < i.end
        ]
        rev = [
            i for i in self.rev_primers
            if i.start < qend and (i.end - 5) <= qend <= (i.end + 5)
        ]
        return tuple(fwd), tuple(rev)

    def primer_check(self, qstart, qend):
        """ Tell if interval fits perfectly inside a primer pair """
        fwd = [i for i in self.fwd_primers if i.end + 1 == qstart]
        rev = [i for i in self.rev_primers if qend == i.start - 1]
        return tuple(fwd), tuple(rev)
