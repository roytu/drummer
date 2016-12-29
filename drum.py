
""" Module for Drum class

    The solution of our drum head is given by a radial field u(r, theta),
    representing the deviation from the "still" drum head.  They take
    the form:

    u_{mn}(r, theta, t) = R(r) * Theta(theta) * T(t)

    where R(r)         = J_m(lamb_{mn} r),
          Theta(theta) = cos (m (theta - theta_0))
          T(t)         = cos (c lamb_{mn} (t - t_0))

    with damping (TODO check this):
    Source: http://www.math.ust.hk/~machas/drum/
          T(t)         = cos (sqrt{ c^2 lamb_{mn}^2 - b^2} (t - t_0))
          T(t)         = e^{-bt} e^{-sqrt{b^2 - lamb_{mn}^2 c^2}t}
    where `b` is the damping factor (units inverse time)

    for m = 0 or positive,
        n = positive
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv, jn_zeros
from scipy.interpolate import griddata
from scipy.io import wavfile

from bessel import Bessel

class Drum(object):
    def __init__(self, a=1, K=0.2, c=10, rcount=30, thetacount=60, m_max=2, n_max=5):
        self.a = a  # Radius
        self.K = K  # Decay factor
        self.c = c  # Tension constant

        self.m_max = m_max
        self.n_max = n_max
        self.rcount = rcount
        self.thetacount = thetacount

        self.t = 0  # Time

        # Generate space
        self.rs = np.linspace(0, a, self.rcount)
        self.thetas = np.linspace(0, 2 * np.pi, self.thetacount)

        self.grid_r, self.grid_theta = np.meshgrid(self.rs, self.thetas)
        self.points = np.array([self.grid_r, self.grid_theta])
        self.points = np.swapaxes(self.points, 0, 1)
        self.points = np.swapaxes(self.points, 1, 2)
        self.points = np.reshape(self.points, (self.rcount * self.thetacount, 2))

        # planes stores the plane wave solutions on impact
        # These should be evolved over time and added together
        self.planes = np.array([])
        self.hit_time = None
        self.ks = np.array([])

    def value(self, ts):
        """ Get the drum state at time `ts` (a numpy array).

            Args:
                ts: times (numpy array of floats) (in seconds)

            Returns:
                numpy array with shape 
        """
        # T(t) = cos (sqrt{ c^2 lamb_{mn}^2 - K^2} (t - t_0)) exp(-Kt)
        # T(t) = cos (c lamb_{mn} (t - t_0))
        if self.planes.size == 0:
            return np.zeros(np.shape(self.points))
        dts = ts - self.hit_time
        
        discriminant = self.K ** 2 - (self.c * self.ks) ** 2 + 0j
        rate = -self.K - np.sqrt(discriminant)
        f1 = np.exp(np.outer(rate, dts)).real

        return np.dot(f1.T, self.planes)

    def hit(self, pos, force):
        """ Simulates a hit of the drum with the given force at some
            position.

            This function takes the impulse to be:

            G = G_0 \delta^{(3)} ( r_0, theta_0, t_0 )

            where t_0 is the current time, and performs a Hankel
            transform to calculate the Bessel coefficients for the
            radial solution.  The theta, time equations are chosen
            such that the function is maximal at the correct
            theta and time.

            This occurs for all modes and the resulting solution
            is added to the current field.

            Force is measured in units of u sec^-2.

            Args:
                pos = (r_0, theta_0, t_0): position and time
                force: impulse strength (in m sec^-2)
        """
        (r_0, theta_0, t_0) = pos

        # Clear arrays
        self.planes = np.array([])
        self.ks = np.array([])

        # Generate a plane
        for m in range(0, self.m_max + 1):
            coeffs = Bessel.coeffs_impulse(m, self.n_max, r_0, self.a, force)

            ks = jn_zeros(m, self.n_max) / self.a
            for k, n, C in zip(ks, range(1, self.n_max + 1), coeffs):
                plane = np.array([])

                rs = self.points[:, 0]
                thetas = self.points[:, 1]

                Rs = jv(m, k * rs)
                Thetas = np.cos(m * (thetas - theta_0))
                plane = C * Rs * Thetas

                # Save the plane
                if self.planes.size == 0:
                    self.planes = np.array([plane])
                else:
                    self.planes = np.vstack([self.planes, plane])

                # Save the hit time and m values
                self.ks = np.append(self.ks, k)
        self.hit_time = t_0

    def wave_value_from_values(self, values, d=10):
        """ Same as wave_value but with the values pre-provided. """
        rs = self.points[:, 0]
        thetas = self.points[:, 1]
        xs = rs * np.cos(thetas)
        ys = rs * np.sin(thetas)
        ds = (d - xs) ** 2 + ys ** 2  # actually distance squared
        amplitudes = np.dot(values, ds) / (d ** 2)

        # TODO normalize?
        #full_amplitude = len(rs) * d  # TODO make this better

        return amplitudes

    def wave_value(self, ts, d=10):
        """ Calculate the amplitude of the wave that an observer hears
            at times `ts`, if the observer is `d` meters away from the
            east side of the drum.

            Args:
                ts: times (numpy array of floats) (in seconds)
                d: distance (in meters)
                
            Returns:
                float
        """
        values = self.value(ts)
        return self.wave_value_from_values(values, d)

    @staticmethod
    def _test_heatmap():
        """ Plots a sample heatmap.
        """
        # Hit the drum at time t = 0
        drum = Drum()
        drum.hit((0.8, 0, 0), 1)

        # Calculate drum values
        t = 0
        points = drum.points
        values = drum.value(0)
        
        # Polar heatmap
        data = griddata(points, values, (drum.grid_r, drum.grid_theta), method='cubic', fill_value=0)
        ax1 = plt.subplot(projection="polar")
        ax1.pcolormesh(drum.thetas, drum.rs, data.T)
        plt.show()

    @staticmethod
    def _test_sound():
        """ Plots a sample sound wave
        """
        # Observer is 10 meters away
        d = 10

        rcount = 30
        thetacount = 60
        m_max = 1
        n_max = 20
        a = 0.5
        c = 150
        K = 10

        # Hit the drum at time t = 0
        drum = Drum(a=a, c=c, K=K, rcount=rcount, thetacount=thetacount, m_max=m_max, n_max=n_max)
        drum.hit((0.8, 0, 0), 1)

        # Calculate drum values
        times = np.linspace(0, 1, 1 * 44100)
        values = drum.wave_value(times, d)

        # Save sound file
        samples = np.array(values / np.max(np.abs(values)), dtype=np.float32)
        wavfile.write("test.wav", 44100, samples)

        # Plot
        fig, ax = plt.subplots()
        ax.plot(times, values)
        plt.show()

if __name__ == "__main__": 
    Drum._test_sound()
