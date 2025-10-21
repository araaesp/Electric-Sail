from numpy import nan

from model.keplerian_body import KeplerianBody


class KeplerianSolarProbe(KeplerianBody):
    def __init__(self, name, mass, a, e, i, omega, OMEGA, radius, alpha, delta, area, M=nan, theta=nan,
                 central_body=None, keplerian_central_body=None,
                 r_diff=0, r_spec=1, e_f=0, e_b=0, a_f=0, a_b=0, chi_f=2 / 3, chi_b=2 / 3):
        """
        Class constructor

        Parameters:
            mass: float. Body mass in kg
            a: semi-major axis in km
            e: eccentricity
            i: inclination (degrees)
            omega: Argument of periapsis (degrees)
            OMEGA: longitude of ascending node (degrees)
            M: mean anomaly (degrees)
            theta: true anomaly (degrees)
            radius: float. Body radius in km
        """
        KeplerianBody.__init__(self, name, mass, a, e, i, omega, OMEGA, radius, M, theta, central_body,
                               keplerian_central_body)
        self.alpha = alpha
        self.delta = delta
        self.area = area
        self.r_diff = r_diff
        self.r_spec = r_spec
        self.e_f = e_f
        self.e_b = e_b
        self.a_f = a_f
        self.a_b = a_b
        self.chi_f = chi_f
        self.chi_b = chi_b
