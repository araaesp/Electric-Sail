import numpy as np
from model.solar_probe import SolarProbe


class SailDynamic:

    @staticmethod
    def calculate_lightness_vector(body):
        """

        :type body: SolarProbe
        """
        sigma_c = 1.5368  # g/m^2
        load_ratio = (sigma_c / (2 * body.sigma))

        lx = load_ratio * ((2 * body.r_spec * np.cos(body.delta) ** 3 * np.cos(body.alpha) ** 3 +
                           (body.chi_f * body.r_diff + body.kappa * body.a_f) * np.cos(body.delta) ** 2 * np.cos(
                    body.alpha) ** 2) + (body.a_f + body.r_diff) * np.cos(body.delta) * np.cos(body.alpha))

        ly = load_ratio * (2 * body.r_spec * np.cos(body.delta) ** 3 * np.cos(body.alpha) ** 2 * np.sin(body.alpha) +
                           (body.chi_f * body.r_diff + body.kappa * body.a_f) * np.cos(body.delta) ** 2 * np.cos(
                    body.alpha) * np.sin(body.alpha))

        lz = load_ratio * (2 * body.r_spec * np.cos(body.alpha) ** 2 * np.cos(body.delta) ** 2 * np.sin(body.delta) +
                           (body.chi_f * body.r_diff + body.kappa * body.a_f) * np.cos(body.delta) * np.cos(
                    body.alpha) * np.sin(body.delta))

        return np.array([lx, ly, lz])

    @staticmethod
    def calculate_acceleration(body):
        """

        :type body: SolarProbe
        """
        # TODO - REMOVE SUN MASS AND G FROM HERE (PUT IN CONSTANT)
        sun_mass = 1.98847e30
        g = 6.6742e-20

        r = body.position
        v = body.velocity
        lightness = SailDynamic.calculate_lightness_vector(body)
        R = np.linalg.norm(r)
        h = np.cross(r, v)
        H = np.linalg.norm(h)

        a_hif_x = (lightness[0] / R * r[0] +
                   lightness[1] / (H * R) * (r[2] * (r[2] * v[0] - r[0] * v[2]) - r[1] * (r[0] * v[1] - r[1] * v[0])) +
                   lightness[2] / H * (r[1] * v[2] - r[2] * v[1]))

        a_hif_y = (lightness[0] / R * r[1] +
                   lightness[1] / (H * R) * (r[0] * (r[0] * v[1] - r[1] * v[0]) - r[2] * (r[1] * v[2] - r[2] * v[1])) +
                   lightness[2] / H * (r[2] * v[0] - r[0] * v[2]))

        a_hif_z = (lightness[0] / R * r[2] +
                   lightness[1] / (H * R) * (r[1] * (r[1] * v[2] - r[2] * v[1]) - r[0] * (r[2] * v[0] - r[0] * v[2])) +
                   lightness[2] / H * (r[0] * v[1] - r[1] * v[0]))

        return sun_mass * g / R ** 2 * np.array([a_hif_x, a_hif_y, a_hif_z])
