import numpy as np
from scipy.integrate import RK45
from scipy.integrate import solve_ivp


class TwoBodySolver:

    def __init__(self, m1, m2):
        # gravitational constant km^3/(kg*s^2)
        self.__g = 6.6742e-20
        self.__m1 = m1
        self.__m2 = m2

    def calculate_derivatives(self, time, state_vector):
        state_vector = np.array(state_vector).reshape(4, 3)
        r1 = state_vector[0]
        r2 = state_vector[1]
        v1 = state_vector[2]
        v2 = state_vector[3]

        r = np.linalg.norm(r2 - r1)

        a1 = self.__g * self.__m2 * (r2 - r1) / r ** 3
        a2 = self.__g * self.__m1 * (r1 - r2) / r ** 3

        return np.array([v1, v2, a1, a2]).flat

    def solve(self, r1, r2, v1, v2, time_span):
        # inputs
        # r1 = initial position vector of body 1
        # r2 = initial position vector of body 2
        # v1 = initial velocity vector of body 1
        # v2 = initial velocity vector of body 2
        #
        # r1 = np.array([0, 3000, 0])
        # r2 = np.array([0, 0, 0])
        #
        # v1 = np.array([70, 0, 0])
        # v2 = np.array([0, 0, 0])

        y0 = np.array([r1, r2, v1, v2]).flat

        def calculate_colision(time, y):
            return abs(y[0] - y[3]) + abs(y[1] - y[4]) + abs(y[2] - y[5])

        calculate_colision.terminal = True

        # integrator = RK45(calculate_derivatives, t0=0.00, y0=y0, t_bound=3600*24/1)
        # integrator.step()

        return solve_ivp(self.calculate_derivatives, (0, time_span), y0, max_step=5.0, events=calculate_colision)

    # output = calculate_derivatives(0.0, y0)
    # print(output)
