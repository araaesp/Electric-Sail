import numpy as np
import os

from solar_sail.sail_dynamic import SailDynamic
from model.body import Body
from model.solar_probe import SolarProbe
from scipy.integrate import RK45
from itertools import permutations
from utils.decorators import log_time
from model.electric_probe import ElectricSailProbe
from solar_sail.electric_sail_dynamic import ElectricSailDynamic


def _check_collision(body1, body2):
    if np.linalg.norm(body1.position - body2.position) < body1.radius + body2.radius:
        raise ArithmeticError('Collision detected between {} and {}'.format(body1.name, body2.name))


class CollisionException(Exception):
    def __init__(self, message, simulation):
        self.message = message
        self.simulation = simulation


class NBodySolver:

    def __init__(self, bodies):
        """
        Class Constructor

        The solver when instantiated will store a list of bodies in its
        private property.
        
        Parameters:
        bodies: list[Body]
            Body is a type declared in this project. It contains several properties
            that will be used during the process. (see body.py)
            
        Class Properties:
        g: private float
            The gravitational constant in km^3/(kg*s^2)
        bodies: private list[Body]
            the bodies that will be calculated in the simulation
        n_bodies: private int
            size of the list of bodies
        """

        self.__g = 6.6742e-20
        self.__bodies = bodies
        self.__n_bodies = len(self.__bodies)
        with open('gravity_log.csv', 'w') as f:
            f.write('time,a_x,a_y,a_z\n')

    @property
    def bodies(self):
        return self.__bodies

    @property
    def n_bodies(self):
        return len(self.bodies)

    def __update_bodies(self, state_vector):
        """
        Private void method to update the properties of each body given a state vector

        Parameters:
        state_vector: array_like, shape(n_bodies + 1,)
            vector with position and velocity of each object
            Ex: [X1, Y1, Z1, ... Xn, Yn, Zn, Vx1, Vy1, Vz1... Vxn, Vyn, Vzn]
        """

        # reshape the vector to be a matrix (n*2 x 3)
        # then update the position and velocity of each body with the current state
        # first half (up/down) in the matrix is position and second half is velocity
        y = state_vector.reshape(self.n_bodies * 2, 3)
        for i in range(self.n_bodies):
            self.bodies[i].position = y[i]
            self.bodies[i].velocity = y[i + self.n_bodies]

    def __calc_g_acceleration(self, body1, body2):
        """
        Private Method
        Calculates Newton gravitational law F = (G*M1*M2/R^3)*r

        Parameters:
        body1: first body object (must have position and mass attributes)
        body2: second body object (must have position and mass attributes)

        return array of accelerations [Ax12 Ay12 Az12]
        """
        # Norm of vector r -> np.lingalg.norm(r)
        r = body2.position - body1.position
        return self.__g * body2.mass * r / np.linalg.norm(r) ** 3

    def calculate_derivatives(self, t, state_vector):
        """
        Connects the scipy RK45 with the calculate_g_acceleration method and solar sail dynamics

        This method delegates the bodies and saves the result of __calc_g_acceleration.
        It takes the number of bodies in the instance, update them with the current state vector
        and calculate the acceleration for each possible peer. Then assembles and returns the 
        time derivative of the state vector.

        Parameters:
        time: current step time. Not used, but it's requested by the RK45
        state_vector: array_like, shape(n_bodies + 1,)
            vector with position and velocity of each object
            Ex: [X1, Y1, Z1, ... Xn, Yn, Zn, Vx1, Vy1, Vz1... Vxn, Vyn, Vzn]

        return time derivative of the state vector.
            [Vx1, Vy1, Vz1... Vxn, Vyn, Vzn, Ax1, Ay1, Az1, ... Axn, Ayn, Azn,]
        """

        n = self.n_bodies
        n_p = n - 1  # number of possible peers for 1 body

        # a_components: 2-D Array with each line being the acceleration [Axmn Aymn Azmn] of each body m caused by a
        # possible peer n. The number of lines for each body will be = n_p
        a_components = np.empty((0, 3))

        # a_sum: array with the sum of the lines of a_components for each body (final acceleration components)
        a_sum = np.empty((0, 3))

        self.__update_bodies(state_vector)

        # for each possible peer of bodies, calculate the acceleration and append it in the a_components array
        for body1, body2 in permutations(self.bodies, 2):
            _check_collision(body1, body2)
            a_components = np.append(a_components, [self.__calc_g_acceleration(body1, body2)], axis=0)

        # for each body, sum the respective lines of a_components and append it in the a vector
        # Ex: for body 2 (n=1) in a 4 body system the components will start in the 4th position and go 3 positions ahead
        # hence start index ->1*3 = 3 and last index -> 1*3 + 3 = 6 (will retrieve elements 3,4,5 of the array).
        for i in range(n): 
            body_n_components = a_components[i * n_p: i * n_p + n_p]
  
            if isinstance(self.bodies[i], SolarProbe):
                a_rad_pressure = SailDynamic.calculate_acceleration(self.bodies[i])
                body_n_components = np.append(body_n_components, [a_rad_pressure], axis=0)

            if isinstance(self.bodies[i], ElectricSailProbe):
            #if i == 1:
                a_electric_sail = ElectricSailDynamic.calculate_acceleration(self.bodies[i])
                body_n_components = np.append(body_n_components, [a_electric_sail], axis=0)

            a_sum = np.append(a_sum, [np.sum(body_n_components, axis=0)], axis=0)

        # to assemble the return, take the velocities from the state vector and append with the accelerations a_sum.
        # the velocities will be after the 3 coordinates of the bodies so they start at position n*3 of the state vector
        return np.append(state_vector[n * 3:], a_sum)

    @log_time
    def solve(self, time_span, t0=0.0, max_step=5.0):
        """
        NBodySolver Callable function

        This method assembles the initial state vector y0 and the integrator, then proceeds to do the integration steps 
        until the time reaches the time_span.

        Parameters:
        time_span: float. period of time (in seconds) that will be ....
        t0: float. initial time (in seconds). Default = 0 s
        max_step: maximum time between each step of the integrator. Default = 5 s

        :returns the array y_t containing the state vectors for each time of the integration. time in seconds and
        coordinates in km
            [
                [t0, X1, Y1, Z1, ... Xn, Yn, Zn, Vx1, Vy1, Vz1... Vxn, Vyn, Vzn],
                [t1, X1, Y1, Z1, ... Xn, Yn, Zn, Vx1, Vy1, Vz1... Vxn, Vyn, Vzn],
                .
                .
                .
                [t_end, X1, Y1, Z1, ... Xn, Yn, Zn, Vx1, Vy1, Vz1... Vxn, Vyn, Vzn],
            ]
        """

        y0 = np.empty((0, 0))  # initial state vector
        t = t0  # setting time to initial time

        # filling initial state vector with each body position
        for body in self.bodies:
            y0 = np.append(y0, body.position)

        # filling initial state vector with each body position
        for body in self.bodies:
            y0 = np.append(y0, body.velocity)

        # initialize empty y_t with number of columns equals to time + 3 coordinates and 3 velocities for each body
        y_t = np.empty((0, 1 + self.n_bodies * 3 * 2))
        # creating the first line of y_t then appending it to y_t
        y_t = np.append(y_t, [np.insert(y0, 0, t0)], axis=0)

        # create the first line of the csv
        with open('gravity_log.csv', 'a') as f:
            f.write('0,0,0,0\n')

        # creating the runge-kutta 4-5 integrator
        integrator = RK45(lambda t, y: self.calculate_derivatives(t, y), t0=t0, y0=y0, t_bound=t0 + time_span, max_step=max_step)

        # doing the steps for the time_span asked and appending each moment to y_t
        while t < t0 + time_span:
            try:
                integrator.step()
            except ArithmeticError as e:
                raise CollisionException(e, y_t)

            # updating the class bodies with current position
            self.__update_bodies(integrator.y)
            # creating current y_t then appending it to y_t
            y_t = np.append(y_t, [np.insert(integrator.y, 0, integrator.t)], axis=0)
            # updating t with the t at the end of the step
            t = integrator.t

            # --- INÍCIO DA NOVA MODIFICAÇÃO ---
            # Após guardar o estado, vamos calcular a gravidade total na nossa sonda (corpo 1)
            # REUTILIZANDO a função de cálculo de gravidade já existente.
            
            sonda = self.__bodies[1]
            a_grav_total_na_sonda = np.array([0.0, 0.0, 0.0])

            # Itera sobre todos os corpos para somar a sua gravidade na sonda
            for outro_corpo in self.__bodies:
                if outro_corpo is not sonda: # Um corpo não exerce gravidade sobre si mesmo
                    # Chama a função interna do solver para calcular a gravidade que 'outro_corpo' exerce na 'sonda'
                    a_grav_total_na_sonda += self.__calc_g_acceleration(sonda, outro_corpo)

            # Guarda o tempo e o vetor de aceleração no ficheiro
            with open('gravity_log.csv', 'a') as f:
                f.write(f'{t},{a_grav_total_na_sonda[0]},{a_grav_total_na_sonda[1]},{a_grav_total_na_sonda[2]}\n')
            # --- FIM DA NOVA MODIFICAÇÃO ---

        return y_t


if __name__ == '__main__':
    print("this file is not for execution")
