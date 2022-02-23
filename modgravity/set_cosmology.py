import numpy as np
import camb


class set_cosmology:
    def __init__(self):
        self.pars = camb.CAMBparams()
        self.pars.set_cosmology(H0=72, ombh2=.0223828, omch2=0.1201075, omk=0, tau=0.06)
        self.results = camb.get_background(self.pars)
        self.z = z

    def get_pars(self):
        return self.pars

    def get_results(self):
        return self.results

    def get_bg(self):
        return self.results.calc_background(self.pars)

    def get_ld(self, z):
        return self.results.luminosity_distance(z)

    def get_hp(self, z):
        return self.results.hubble_parameter(z)
