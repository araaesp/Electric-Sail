import numpy as np
from datetime import datetime
from pathlib import Path
from model.electric_probe import ElectricSailProbe

K_T = 4.3          # Coeficiente do Hoytether (ex: para 4 sub-fios)
BETA = 0.25        # Razão massa-potência da espaçonave (kg/W)
N_EARTH = 7.3e6    # Densidade de elétrons do vento solar a 1 UA (partículas/m^3)
RHO_W = 4000       # Massa específica do material do fio (kg/m^3)
E_CARGA = 1.602e-19  # Carga elementar (C)
M_ELETRON = 9.109e-31 # Massa do elétron (kg)


class ElectricSailDynamic:
    _v_sw_km_s = 400.0
    _n_earth_m3 = 7.3e6
    _t_e_earth_ev = 10.0
    _config_source = "default"

    @classmethod
    def configure_from_csv(cls, start_time, csv_path=None):
        if csv_path is None:
            csv_path = Path(__file__).resolve().parent.parent / "atividade_solar.csv"

        data = np.loadtxt(csv_path)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        if data.shape[1] < 6:
            raise ValueError("atividade_solar.csv must have 6 columns")

        if isinstance(start_time, str):
            start_time = datetime.strptime(start_time, "%Y-%m-%d %H:%M:%S")

        year = int(start_time.year)
        doy = int(start_time.timetuple().tm_yday)
        hour = int(start_time.hour)

        mask_day = (data[:, 0].astype(int) == year) & (data[:, 1].astype(int) == doy)
        if not np.any(mask_day):
            raise ValueError(f"No solar wind row found for year={year}, doy={doy}")

        candidates = data[mask_day]
        hour_col = candidates[:, 2]
        idx = int(np.argmin(np.abs(hour_col - hour)))
        row = candidates[idx]

        temp_k = float(row[3])
        dens_cm3 = float(row[4])
        vel_km_s = float(row[5])

        cls._v_sw_km_s = vel_km_s
        cls._n_earth_m3 = dens_cm3 * 1e6
        cls._t_e_earth_ev = temp_k / 11604.5
        cls._config_source = f"{year}-DOY{doy}-H{int(round(row[2]))}"

    @staticmethod
    def calculate_thrust_per_m(body, r_m: float = None):
        """
        Calcula a força de empuxo por metro de fio (σ_F) em N/m.

        Se r_m (distância ao Sol em metros) for dada, n e T_e são
        calculados em funçao de r:
            n   = n_earth * (r_earth / r)^2          (eq 68)
            T_e = T_e_earth * (r_earth / r)^(1/3)    (eq 69)
        Se não tiver r_m, usa os valores constantes a 1 U que ja tem.
        """
        m_p = 1.6726219e-27  # Massa do proton em kg
        epsilon_0 = 8.854187817e-12 # Permissividade do vacuo
        e = 1.60217662e-19

        r_earth = 1.496e11   # 1 UA em m

        v_sw = ElectricSailDynamic._v_sw_km_s * 1e3
        n_earth = ElectricSailDynamic._n_earth_m3
        T_e_earth_eV = ElectricSailDynamic._t_e_earth_ev

        if r_m is not None and r_m > 0:
            # eq. 68
            n = n_earth * (r_earth / r_m) ** 2
            # eq. 69
            T_e_joules = (T_e_earth_eV * (r_earth / r_m) ** (1.0 / 3.0)) * e
        else:
            # valores a 1 UA
            n = n_earth
            T_e_joules = T_e_earth_eV * e

        # Parametros da vela
        V = body.V
        r_w = body.r_w * 1000

        termo_raiz_ln = np.sqrt((epsilon_0 * T_e_joules) / (n * e**2))
        arg_ln = (2 / r_w) * termo_raiz_ln
        resultado_ln = np.log(arg_ln)
        expoente = (m_p * v_sw**2 / (e * V)) * resultado_ln
        denominador = e * np.sqrt(np.exp(expoente) - 1)
        numerador = 6.18 * m_p * v_sw**2 * np.sqrt(n * epsilon_0 * T_e_joules)

        sigma_F = numerador / denominador

        return sigma_F

    @staticmethod
    def calculate_acceleration(body: ElectricSailProbe):
        """
        Calcula o vetor de aceleração gerado pela Vela Elétrica.
        """
        r_base_UA = 1.496e11

        # Parâmetros da vela
        
        N = body.N
        L = body.L * 1000
        phi = body.phi
        theta = body.theta
        r_w = body.r_w * 1000
        V = body.V

        # Posicao no referencial inercial
        r_inercial = body.position * 1000
        v_inercial = body.velocity * 1000
        
        # Distancia ao Sol
        r = np.linalg.norm(r_inercial)

        # Empuxo/m calculado na distância atual (n e T_e variam com r)
        sigma_F_base = ElectricSailDynamic.calculate_thrust_per_m(body, r_m=r)

        magnitude_forca = (1/2) * N * L * sigma_F_base * (r_base_UA / r)**(7/6)
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)

        Fx_orbita = cos_phi * sin_theta * cos_theta
        Fy_orbita = -sin_phi * cos_phi * cos_theta**2
        Fz_orbita = cos_phi**2 * cos_theta**2 + 1
        
        F_vela_orbita = magnitude_forca * np.array([Fx_orbita, Fy_orbita, Fz_orbita])
        #F_vela_orbita = 1e-6 * np.array([1, 0, 0])

        # eixo z da orbita escrito no sistema inercial
        z_o = r_inercial / r

        # eixo x da orbita escrito no sistema inercial
        h_vec = np.cross(r_inercial, v_inercial)
        x_o = - h_vec / np.linalg.norm(h_vec)

        # eixo y da orbita escrito no sistema inercial
        y_o = np.cross(z_o, x_o)

        # Matriz de rotacao orbita -> inercial
        matriz_rotacao = np.array([x_o, y_o, z_o]).T

        # rotacao
        F_vela_inercial = matriz_rotacao @ F_vela_orbita

        # Massa do corpo da vela por metro do fio (eq 75)
        sigma_mb = 2 * K_T * BETA * N_EARTH * r_w * np.sqrt((2 * E_CARGA**3 * V**3) / M_ELETRON)
        
        # Massa dos fios por metro de fio (eq 76)
        sigma_mt = K_T * np.pi * RHO_W * r_w**2
        
        # Massa seca total por metro de fio
        sigma_m_seca = sigma_mb + sigma_mt
        
        L_total = N * L
        m_seca = sigma_m_seca * L_total
        
        # a = F/m
        a_vela_inercial = (F_vela_inercial / m_seca) * 0.001  # converter de m/s^2 para km/s^2

        return a_vela_inercial