from model.body import Body


class SolarProbe(Body):
    def __init__(self, name, mass, position, velocity, radius,
                 alpha, delta, area,
                 r_diff=0, r_spec=1, e_f=0, e_b=0, a_f=0, a_b=0, chi_f=2/3, chi_b=2/3):
        Body.__init__(self, name, mass, position, velocity, radius)
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

    @property
    def sigma(self):
        return self.mass*1000/self.area     # g/m^2

    @property
    def kappa(self):
        if self.e_b == 0 and self.e_f == 0:
            return 0
        else:
            return (self.chi_f * self.e_f - self.chi_b * self.e_b) / (self.e_f + self.e_b)
