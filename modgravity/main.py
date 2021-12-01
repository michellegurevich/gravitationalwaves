import numpy as np
import math
from matplotlib import pyplot as plt

from SetCosmology import SetCosmology
from CalculateDistances import CalculateDistances
from ModifiedPolarization import ModifiedPolarization
from TimeDomain import TimeDomain
from Plots import Plots
from Model import Model


def main():
    TD = TimeDomain()
    MP = ModifiedPolarization()
    M = Model()

    TaylorF2 = {
        'approximant': 'TaylorF2',
        'mass1': 10, 'mass2': 15,
        'delta_f': 1.0 / 170, 'f_lower': 40,
    }

    TaylorT2 = {
        'approximant': 'TaylorT2',
        'mass1': 10, 'mass2': 15,
        'delta_t': 1.0 / 4096,
        'f_lower': 40
    }

    IMRPhenomA = {
        'approximant': 'IMRPhenomA',
        'mass1': 65, 'mass2': 80,
        'delta_f': 1.0 / 4, 'f_lower': 40,
        'delta_t': 1.0 / 4096
    }

    test = {
        'approximant': 'test',
        'mass1': TD.m_1, 'mass2': TD.m_2,
        'delta_f': 1.0 / 320, 'f_lower': 40,
        'delta_t': 1.0 / 4096
    }

    wf = M.get_waveform()
    f, z_max, z = M.get_params(wf)
    #print(wf)
    #print(f, z_max, z)
    amplitude, std_phase = M.decompose_waveform(wf, f, z_max, z)
    print(amplitude)
    print(std_phase)

    alpha = 4
    A_term = .001

    mod_amplitude, mod_phase = M.perform_modification(wf, std_phase, f, z_max, z, alpha, A_term)
    print(mod_amplitude, mod_phase)

if __name__ == '__main__':
    main()
