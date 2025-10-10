from enum import IntEnum

from .hmm import HMM
from .hmm_model_data import HMM_MODEL_DATA


models = {
    name: HMM(name, **data)
    for name, data
    in HMM_MODEL_DATA.items()
}
""" the HMM models """


Err = IntEnum('Err', ' '.join((f'E{i}' for i in range(20))))
""" common error codes """
