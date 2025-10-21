from enum import IntEnum

from .hmm import HMM
from .hmm_model_data import HMM_MODEL_DATA
from .primers import Primer


Err = IntEnum('Err', ' '.join((f'E{i}' for i in range(20))))
""" common error codes """


def get_models():
    """ Get known HMM models """
    primers = Primer.load()
    return {
        name:
        HMM(name, **data, primers=[i for i in primers if i.hmm_name == name])
        for name, data
        in HMM_MODEL_DATA.items()
    }
