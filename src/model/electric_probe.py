from model.body import Body
import numpy as np

class ElectricSailProbe(Body):
    def __init__(self, name, mass, position, velocity, radius,
                 N, L, V, r_w, phi_deg, theta_deg):
        super().__init__(name, mass, position, velocity, radius)
        self.N = N              # Número de fios
        self.L = L              # Comprimento de cada fio (m)
        self.V = V              # Tensão da vela (V)
        self.r_w = r_w          # Raio do fio (m)
        self.phi = np.deg2rad(phi_deg)    # Ângulo phi (radianos)
        self.theta = np.deg2rad(theta_deg) # Ângulo theta (radianos)