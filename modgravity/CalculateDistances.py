import numpy as np

from SetCosmology import SetCosmology


class CalculateDistances:
    SC = SetCosmology()

    def set_redshift(self, z):
        pass

    def get_redshift(self):
        return self.z

    def alpha_dist(self, z_max, a):
        return pass

    def lum_dist(self, z_max):
        # FROM eqn (2) THERE IS FACTOR OF POWERS OF C IN ENERGY TERM ->

        def func(z_prime):
            return 1 / np.sqrt(z_prime)  # issue because math sqrt fctn converts array into scalar, use np.sqrt

        def func2(z_prime):
            return omega_m * (1 + a * z_prime) ** 3 + omega_lambda

        def integrand(z_prime):
            return func(func2(z_prime))

        a = np.array(range(0, 50))  # controls number of points integrand is evaluated at
        integral = simps(integrand(z_prime), z_prime, axis=0)

        def lum_dist(z):
            # constant_term = c * (1 + z) / h0
            constant_term = (1 + z) / h0
            # return constant_term * integral
            return np.flip(constant_term * integral)  # graviton travels from source to observer#
