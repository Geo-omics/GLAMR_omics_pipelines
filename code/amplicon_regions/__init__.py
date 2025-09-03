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
    # mito 16S:
    # from set_15:
    #       16S AR start:848
    #       16S BR end:1474
    ('mitochondria', '16S_rRNA'): {
        'AR_BR': (840, 1480),
    },
    ('archaea', '16S_rRNA'): ('bacteria', '16S_rRNA'),  # TODO
    ('eukaryote', '18S_rRNA'): {
        # adapted from glamr_data_management repository
        # per lit. F04/R22 cover V1-V2
        # FIXME
        'V1': (40, 199),  # SSU_F04 starts at 30, end uncertain
        'V2': (200, 420),  # start uncertain, SSU_R22 ends at 435
        'V3': (500, 699),  # start < 600
        'V4': (700, 899),  # start < 800
        'V5': (1000, 1199),  # start < 1100
        'V6': (1300, 1399),  # start < 1350
        'V7': (1400, 1599),  # start < 1500
        'V8': (1600, 1699),  # start < 1650
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
