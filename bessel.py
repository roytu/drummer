
""" Module for Bessel function calculations """

import numpy as np
from scipy.special import jv, jn_zeros
import matplotlib.pyplot as plt

class Bessel(object):
    @staticmethod
    def coeffs_impulse(m, ncount, r_0, a, G):
        """ Calculate the coefficients for the Bessel function `m`, `n`
            with an impulse at `r_0` and a radius of `a`.

            We should be able to recover the original function with:

                f(x) = sum_{n=1}^inf c_n J_m((alpha_{mn} / a) x)

            where alpha_{mn} is the root of the Bessel function.

            Source: http://www.hit.ac.il/staff/benzionS/Differential.Equations/Orthogonality_of_Bessel_functions.htm

            Args:
                m: order of Bessel function
                ncount: number of coeffs to compute
                r_0: radius of impulse
                a: radius of the boundary
                G: strength of impulse
        """
        m = np.float(m)
        a = np.float(a)

        alphas = jn_zeros(m, ncount)

        numer = G * jv(m, alphas / a * r_0) * r_0
        denom = (a ** 2 / 2) * (jv(m + 1, alphas)) ** 2

        return numer / denom

    def coeffs_func(m, coeffs):
        """ TODO doc """
        pass

if __name__ == "__main__":
    m = 0
    ncount = 50
    r_0 = 0.8
    a = 1.0
    steps = 1000

    # Build space and coefficients
    rs = np.linspace(0, a, steps)
    coeffs = Bessel.coeffs_impulse(m, ncount, r_0, a, 1)

    # Evaluate Bessel functions
    bessels = []

    ks = jn_zeros(m, ncount) / a
    for k, n in zip(ks, xrange(1, ncount + 1)):
        ys_ = jv(m, k * rs)
        bessels.append(ys_)

    bessels = np.array(bessels)

    # Add Bessels linearly
    ys = np.dot(coeffs, bessels)

    # Plot should be a Dirac spike at r = r_0
    fig, ax = plt.subplots()
    ax.plot(rs, ys)
    plt.show()
