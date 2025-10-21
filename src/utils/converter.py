import json

import numpy as np
from model.body import Body
from model.keplerian_body import KeplerianBody
from scipy import optimize

from model.solar_probe import SolarProbe
from utils.decorators import log_time

JSON_PATH = '../resources/solar-system.json'


def rotation_z(alpha):
    return np.array([[np.cos(alpha), np.sin(alpha), 0],
                     [-np.sin(alpha), np.cos(alpha), 0],
                     [0, 0, 1]])


def rotation_x(alpha):
    return np.array([[1, 0, 0],
                     [0, np.cos(alpha), np.sin(alpha)],
                     [0, -np.sin(alpha), np.cos(alpha)]])


def get_rotation_matrix(keplerian_body):
    omega_rad = np.deg2rad(keplerian_body.omega)
    i_rad = np.deg2rad(keplerian_body.i)
    OMEGA_RAD = np.deg2rad(keplerian_body.OMEGA)

    r1 = np.dot(rotation_z(-OMEGA_RAD), rotation_x(-i_rad))
    r2 = np.dot(r1, rotation_z(-omega_rad))
    return r2


@log_time
def kepler_to_cartesian(keplerian_body):
    """

    :param keplerian_body: KeplerianBody
    :return: Body
    """
    if keplerian_body.central_body is None:
        keplerian_body.central_body = kepler_to_cartesian(keplerian_body.keplerian_central_body)

    with open(JSON_PATH) as json_file:
        data = json.load(json_file)
        astra = data[keplerian_body.central_body.name.lower()]

    mu = astra['mu']

    # Mean motion
    n = np.sqrt(mu / keplerian_body.a ** 3)

    if np.isnan(keplerian_body.M):
        theta = np.deg2rad(keplerian_body.theta)

        U = np.arctan2(np.sqrt(1 - keplerian_body.e ** 2) * np.sin(theta), keplerian_body.e + np.cos(theta))
    else:
        # converting main anomaly to radians
        M = np.deg2rad(keplerian_body.M)

        # Eccentric anomaly U
        U_function = lambda u: u - keplerian_body.e * np.sin(u) - M
        U_prime = lambda u: 1 - keplerian_body.e * np.cos(u)

        # applying newton-raphson to find U
        U = optimize.newton(func=U_function, x0=M, fprime=U_prime, tol=1e-8)

        # True anomaly theta
        theta = 2 * np.arctan2(np.sqrt(1 - keplerian_body.e) * np.sin(U / 2),
                               np.sqrt(1 - keplerian_body.e) * np.cos(U / 2))

    # Get central body distance
    rc = keplerian_body.a * (1 - keplerian_body.e * np.cos(U))

    # Position and velocity on local coordinates system
    position = [keplerian_body.a * (np.cos(U) - keplerian_body.e),
                keplerian_body.a * np.sin(U) * np.sqrt(1 - keplerian_body.e ** 2),
                0]

    velocity = np.array([-n * keplerian_body.a ** 2 / rc * np.sin(U),
                         n * keplerian_body.a ** 2 / rc * np.cos(U) * np.sqrt(1 - keplerian_body.e ** 2),
                         0])

    # Rotate the coordinate system to the inertial frame.
    # Rotation matrix = R(-OMEGA)*R(-i)*R(-omega)
    rotation_matrix = get_rotation_matrix(keplerian_body)

    # Position and velocity in relation to the central body of the 2-body system
    relative_position = np.dot(rotation_matrix, position)
    relative_velocity = np.dot(rotation_matrix, velocity)

    # Position and velocity in relation to the inertial axis of the whole system.
    absolute_position = relative_position + keplerian_body.central_body.position
    absolute_velocity = relative_velocity + keplerian_body.central_body.velocity

    return Body(keplerian_body.name, keplerian_body.mass, absolute_position, absolute_velocity, keplerian_body.radius)


def probe_kepler_to_cartesian(keplerian_solar_probe):
    """
    :param keplerian_solar_probe: KeplerianSolarProbe
    :return: SolarProbe
    """

    body = kepler_to_cartesian(keplerian_solar_probe)

    return SolarProbe(
        body.name, body.mass, body.position, body.velocity, body.radius,
        keplerian_solar_probe.alpha,
        keplerian_solar_probe.delta,
        keplerian_solar_probe.area,
        keplerian_solar_probe.r_diff,
        keplerian_solar_probe.r_spec,
        keplerian_solar_probe.e_f,
        keplerian_solar_probe.e_b,
        keplerian_solar_probe.a_f,
        keplerian_solar_probe.a_b,
        keplerian_solar_probe.chi_f,
        keplerian_solar_probe.chi_b, )


if __name__ == '__main__':
    central_body = Body('earth', 5.97237e24, [0, 0, 0], [0, 0, 0], 6378.1366)

    example = KeplerianBody('HUBBLE', 10, 6.922453751349309e3, 1.1591437e-3, 2.858871232061450E+01,
                            4.804080043054086E+01,
                            8.037273861036857E+01, 42, central_body, theta=3.321400340149752E+02)
    print(example, kepler_to_cartesian(example), sep='\n===================\n')


def cartesian_to_keplerian(r_vec, v_vec, mu):
    """
    Converte vetores cartesianos (posição e velocidade)
    para elementos orbitais keplerianos. Ajustada para quando a orbita é equatorial ou circular
    """
    r_mag = np.linalg.norm(r_vec)
    if r_mag == 0: return (0,0,0,0,0,0)

    v_mag = np.linalg.norm(v_vec)
    h_vec = np.cross(r_vec, v_vec)
    h_mag = np.linalg.norm(h_vec)
    k_hat = np.array([0, 0, 1])
    n_vec = np.cross(k_hat, h_vec)
    n_mag = np.linalg.norm(n_vec)
    e_vec = (1/mu) * ((v_mag**2 - mu/r_mag) * r_vec - np.dot(r_vec, v_vec) * v_vec)
    e = np.linalg.norm(e_vec)
    
    
    # energia e semieixo maior
    energia = (v_mag**2 / 2) - (mu / r_mag)
    if not np.isclose(e, 1.0):
        a = -mu / (2 * energia)
    else:
        a = np.inf

    # inclinação
    i = np.arccos(h_vec[2] / h_mag)

    # i ≈ 0
    if np.isclose(i, 0.0):
        OMEGA = 0 
        omega = np.arccos(e_vec[0] / e)
        if e_vec[1] < 0:
            omega = 2 * np.pi - omega
    # i ≠ 0
    else:
        OMEGA = np.arccos(n_vec[0] / n_mag)
        if n_vec[1] < 0:
            OMEGA = 2 * np.pi - OMEGA
        
        omega = np.arccos(np.dot(n_vec, e_vec) / (n_mag * e))
        if e_vec[2] < 0:
            omega = 2 * np.pi - omega

    # e ≈ 0
    if np.isclose(e, 0.0):
        if np.isclose(i, 0.0):
            nu = np.arccos(r_vec[0] / r_mag)
            if r_vec[1] < 0:
                nu = 2 * np.pi - nu
        else:
            nu = np.arccos(np.dot(n_vec, r_vec) / (n_mag * r_mag))
            if r_vec[2] < 0:
                nu = 2 * np.pi - nu
    # e ≠ 0 e i ≠ 0
    else:
        nu = np.arccos(np.dot(e_vec, r_vec) / (e * r_mag))
        if np.dot(r_vec, v_vec) < 0:
            nu = 2 * np.pi - nu

    i_deg = np.rad2deg(i)
    OMEGA_deg = np.rad2deg(OMEGA)
    omega_deg = np.rad2deg(omega)
    nu_deg = np.rad2deg(nu)

    return (a, e, i_deg, OMEGA_deg, omega_deg, nu_deg)