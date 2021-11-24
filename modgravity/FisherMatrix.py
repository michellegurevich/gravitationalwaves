import camb
import numpy as np
from SetCosmology import SetCosmology
from ModifiedPolarization import ModifiedPolarization
from TimeDomain import TimeDomain


class FisherMatrix:
    SC = SetCosmology()
    MP = ModifiedPolarization()
    TD = TimeDomain()

    def __init__(self):
        pass

    @classmethod
    def generate_matrix(cls):
        return 0

    @classmethod
    def plot(cls):
        return 0