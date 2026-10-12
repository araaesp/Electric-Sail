from model.body import Body
import numpy as np

class ElectricSailProbe(Body):
    """Sonda e-sail.

    A entrada usa L e r_w em km, V em volts e phi_deg/theta_deg em graus.
    Os ângulos são convertidos para radianos e armazenados como phi/theta.
    """

    def __init__(self, name, mass, position, velocity, radius,
                 N, L, V, r_w, phi_deg, theta_deg):
        super().__init__(name, mass, position, velocity, radius)
        self.N = N              # Número de fios
        self.L = L              # Comprimento de cada fio (km)
        self.V = V              # Tensão elétrica da vela (V)
        self.r_w = r_w          # Raio do fio (km)
        self.phi = np.deg2rad(phi_deg)    # Entrada em graus; armazenado em radianos
        self.theta = np.deg2rad(theta_deg) # Entrada em graus; armazenado em radianos
