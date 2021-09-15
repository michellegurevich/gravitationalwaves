from matplotlib import pyplot as plt
import numpy as np

from CalculateDistances import CalculateDistances
from SetCosmology import SetCosmology


class Plots:

    # plot scale factor against redshift
    def scale_factor(self, z):
        SC = SetCosmology()
        hp = SC.get_hp(z)
        plt.plot(z, hp)
        plt.xlabel('$z$')
        plt.ylabel('$H(t)$')
        plt.title('$H(t)$ as a function of redshift')
        return plt.show()

    # plot ratios of (modified) alpha and (standard) luminosity distances against redshift
    def alpha_lum_ratios(self, z_max):
        results = []
        CD = CalculateDistances()

        alpha_values = np.arange(0, 4.5, .5)  # alpha values end at 4, per GR tests paper specs
        for a in alpha_values:
            alpha_lum_ratio = (CD.alpha_dist(z_max, a) / CD.lum_dist(z_max))[1]
            results.append(alpha_lum_ratio)

        plt.plot(np.hstack(alpha_values), np.hstack(results))
        plt.xlabel(r'$\alpha$')
        plt.ylabel(r'$D_{\alpha} / D_L {Mpc}$')
        plt.title('Ratio of alpha to standard luminosity distances')
        return plt.show()
