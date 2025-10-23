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

        self.fwd_primers = []
        self.rev_primers = []
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
                    self.fwd_primers.append(i)
                elif i.is_reverse():
                    self.rev_primers.append(i)
                else:
                    # skip
                    continue

    def __str__(self):
        return self.name

    def __repr__(self):
        return f'<{self.__class__.__name__}: {self.name}>'

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
            i for i in self.rev_primers.items()
            if i.start < qend and (i.end - 5) <= qend <= (i.end + 5)
        ]
        return tuple(fwd), tuple(rev)

    def primer_check(self, qstart, qend):
        """ Tell if interval fits perfectly inside a primer pair """
        fwd = [i for i in self.fwd_primers if i.end + 1 == qstart]
        rev = [i for i in self.rev_primers if qend == i.start - 1]
        return tuple(fwd), tuple(rev)
