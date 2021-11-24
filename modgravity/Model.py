import camb
import numpy as np
from SetCosmology import SetCosmology
from ModifiedPolarization import ModifiedPolarization
from TimeDomain import TimeDomain


class Model:
    SC = SetCosmology()
    MP = ModifiedPolarization()
    TD = TimeDomain()

    def __init__(self):
        pass

    @classmethod
    def get_waveforms(cls):
        return 0

    @classmethod
    def perform_modification(cls):
        return 0
