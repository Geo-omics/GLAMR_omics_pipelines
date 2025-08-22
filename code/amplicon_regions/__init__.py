REGIONS = {
    ('bacteria', '16S_rRNA'): {
        # result from find_16S_regions_in_hmm()
        'V1': (69, 93),
        'V2': (132, 240),
        'V3': (431, 488),
        'V4': (567, 673),
        'V5': (813, 870),
        'V6': (977, 1034),
        'V7': (1108, 1165),
        'V8': (1235, 1286),
        'V9': (1427, 1457),
    },
}


# consistency check at module load time
for key, regions in REGIONS.items():
    if isinstance(regions, dict):
        cur = 0
        for name, (start, end) in regions.items():
            if end <= start:
                raise ValueError('start must be smaller than end')
            if start <= cur:
                raise ValueError(
                    'regions must be listed in increasing coordinates and not '
                    'overlap'
                )
            cur = end
    else:
        if regions not in REGIONS:
            raise ValueError(f'bad reference {regions}, it\'s not declared')


coordinates = {}
for (tax_group, gene), val in REGIONS.items():
    if isinstance(val, tuple):
        # refer to a previously declared regions listing
        try:
            val = REGIONS[val]
        except KeyError as e:
            raise KeyError(
                f'key {val} must be declared in REGIONS'
            ) from e
    coordinates[(tax_group, gene)] = val
