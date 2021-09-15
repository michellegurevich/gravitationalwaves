import numpy as np
import camb


class SetCosmology:
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=72, ombh2=.0223828, omch2=0.1201075, omk=0, tau=0.06)
    results = camb.get_background(pars)
    bg = results.calc_background(pars)
    ld = results.luminosity_distance(z=np.linspace(0, 4))
    hp = results.hubble_parameter(z=np.linspace(0, 4))

    @staticmethod
    def get_pars(self):
        return self.pars

    @staticmethod
    def get_results(self):
        return self.results

    @staticmethod
    def get_bg(self):
        return self.bg

    def get_ld(self, z):
        return self.ld

    def get_hp(self, z):
        return self.hp
