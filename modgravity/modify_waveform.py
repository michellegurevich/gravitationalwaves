import math
import cmath
import numpy as np
from matplotlib import pyplot as plt

from distances import distances

H = 1  # Planck constant, natural units
F_E = 10 ** 4  #Hz
LAMBDA_G = 1.6 * 10 ** 16  # m
EPSILON_LIGO = 1  # constant term for LIGO
# EPSILON_LISA = math.sqrt(3) / 2  # for LISA


class modify_waveform:

    def __init__(self, cosmo_params, phenom_params, wf_params):
        self.cosmo_params = cosmo_params
        self.phenom_params = phenom_params
        self.wf_params = wf_params
        self.dist = distances(cosmo_params, phenom_params)

    def __unpack_params(self, **params):
        A = phenom_params.get('A') if 'A' in phenom_params else print("A not found")
        return params

    def lambda_A(self, **phenom_params):
        self.__unpack_params(**phenom_params)
        lambda_A = H * phenom_params['A'] ** (1 / (phenom_params['alpha'] - 2))
        return lambda_A

    def u(self, **wf_params):
        self.__unpack_params(**wf_params)
        return math.pi * chirp_mass * f

    def beta(self, **cosmo_params):
        self.__unpack_params(**cosmo_params)
        alpha = 0
        # pass through alpha = 0 for D_0 calc
        D_0 = self.dist.mod_luminosity(self.cosmo_params, self.phenom_params)
        return (math.pi ** 2 * D_0 * chirp_mass) / (LAMBDA_G ** 2 * (1 + z[-1]))

    def zeta():
        return 0

    def chirp_mass():
        return 0

    def chirp_mass_e():
        return 0

    def m():
        return 0

    def v_slo():
        """ proportional to velocity at last stable orbit in the strongly relativistic regime (here we use the
        Schwarzschild metric termination point) """
        return 1 / math.sqrt(6)

    def gr_amplitude():
        """ curved A term; has frequency dependence """
        return 0

    def modified_amplitude():
        """ A~(f) term """
        return 0

    def gr_phase():

        return 0

    def gr_phase_numerical_coefficients():
        return 0

    def phase_perturbation():
        return 0

    def modified_phase():
        return 0

    def standard_polarization():
        return 0
