class HMM:
    """
    HMM model
    """
    def __init__(self, name, regions=None, primers=None, length=None):
        if not name:
            raise ValueError('name must be supplied')

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

        if primers is None:
            primers = dict(fwd={}, rev={})
        # check primer coords consistency
        for primer_coords in primers.values():
            for start, end in primer_coords.values():
                if start <= 0 or end <= 0:
                    raise ValueError('primer coords must be positive')
                if end <= start:
                    raise ValueError('start must be smaller than end')
                if length is not None and length < end:
                    raise ValueError('prim end must not be larger than length')

        self.fwd_primers = primers.get('fwd')
        self.rev_primers = primers.get('rev')
        self.length = length

    def __str__(self):
        return self.name

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
        fwd = []
        rev = []
        for name, (start, end) in self.fwd_primers.items():
            if (start - 5) <= qstart <= (start + 5) and qstart < end:
                fwd.append(name)
        for name, (start, end) in self.rev_primers.items():
            if start < qend and (end - 5) <= qend <= (end + 5):
                rev.append(name)
        return tuple(fwd), tuple(rev)

    def primer_check(self, qstart, qend):
        """ Tell if interval fits perfectly inside a primer pair """
        fwd = []
        rev = []
        for name, (start, end) in self.fwd_primers.items():
            if end + 1 == qstart:
                fwd.append(name)
        for name, (start, end) in self.rev_primers.items():
            if qend == start - 1:
                rev.append(name)
        return tuple(fwd), tuple(rev)
