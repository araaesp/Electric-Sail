from numpy import array, nan


class KeplerianBody:
    def __init__(self, name, mass, a, e, i, omega, OMEGA, radius, M=nan, theta=nan, central_body=None,
                 keplerian_central_body=None):
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
        self.name = name
        self.mass = mass
        self.a = a
        self.e = e
        self.i = i
        self.omega = omega
        self.OMEGA = OMEGA
        self.radius = radius
        self.central_body = central_body
        self.keplerian_central_body = keplerian_central_body
        self.M = M
        self.theta = theta

    def __str__(self):
        return "{}\nmass: {}\na: {}\ne: {}\ni: {}\nomega: {}\nOMEGA: {}\nM: {}\nradius:{}\ntheta {}" \
            .format(self.name, self.mass, self.a, self.e, self.i, self.omega, self.OMEGA, self.M, self.radius,
                    self.theta)
