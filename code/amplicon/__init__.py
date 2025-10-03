from .hmm import HMM
from .hmm_model_data import HMM_MODEL_DATA


models = {
    name: HMM(name, **data)
    for name, data
    in HMM_MODEL_DATA.items()
}
