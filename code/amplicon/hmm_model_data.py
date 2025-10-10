""" This file contains a description of all the HMM models """

HMM_MODEL_DATA = {
    '16S_rRNA_bac': {
        'regions': {
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
        'primers': {
            'fwd': {
                '515F': (505, 523),
                '357F': (338, 354),
                # '341F': (?, ?)
            },
            'rev': {
                '806R': (777, 796),
                '926R': (899, 916),
            },
        },
        'length': 1533,
    },
    # mito 16S:
    # from set_15:
    #       16S AR start:848
    #       16S BR end:1474
    '16S_rRNA_mito': {
        'regions': {
            'AR_BR': (840, 1480),
        },
        'primers': {
            'fwd': {
                '16Sar': (848, 867),
            },
            'rev': {
                '16Sbr-H': (1453, 1474),
            },
        },
    },
    '16S_rRNA_arc': {
        'length': 1477,
    },
    '18S_rRNA_euk': {
        'regions': {
            # adapted from glamr_data_management repository
            # per lit. F04/R22 cover V1-V2
            # FIXME
            'V1': (51, 199),  # SSU_F04 starts at 30, end uncertain
            'V2': (200, 416),  # start uncertain, SSU_R22 ends at 435
            'V3': (500, 699),  # start < 600
            'V4': (700, 899),  # start < 800
            'V5': (1000, 1199),  # start < 1100
            'V6': (1300, 1399),  # start < 1350
            'V7': (1400, 1599),  # start < 1500
            'V8': (1600, 1699),  # start < 1650
        },
        'primers': {
            'fwd': {
                'SSU_F04': (30, 50),
            },
            'rev': {
                'SSU_R22': (417, 435),
            },
        },
        'length': 1851,
    },
    '12S_rRNA_mito': {
        'length': 947,
    },
    '23S_rRNA_bac': {'length': 2906},
    '23S_rRNA_arc': {'length': 2978},
    '28S_rRNA_euk': {'length': 3729},
    '5S_rRNA_arc': {},
    '5_8S_rRNA_arc': {},
    '5S_rRNA_bac': {},
    '5S_rRNA_euk': {},
    '5_8S_rRNA_euk': {},
}
