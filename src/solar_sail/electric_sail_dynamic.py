import numpy as np
from datetime import datetime
from pathlib import Path
from model.electric_probe import ElectricSailProbe

# Constantes usadas nas equações do modelo; valores em SI, salvo indicação.
K_T = 4.3          # Coeficiente do Hoytether
BETA = 0.25        # Razão massa-potência da espaçonave (kg/W)
N_EARTH = 7.3e6    # Densidade fixa da massa (eq. 75), distinta da amostra usada no empuxo.
RHO_W = 4000       # Massa específica do material do fio (kg/m^3)
E_CARGA = 1.602e-19  # Valor arredondado usado no cálculo da massa (C)
M_ELETRON = 9.109e-31 # Massa do elétron (kg)
M_PROTON = 1.6726219e-27  # Massa do próton (kg)
EPSILON_0 = 8.854187817e-12  # Permissividade do vácuo (F/m)
E_CARGA_F = 1.60217662e-19  # Valor mais preciso usado no cálculo da força (C)
AU_M = 1.496e11     # UA (m)


class ElectricSailDynamic:
    """Modelo de empuxo da e-sail e sua transformação para o referencial inercial."""

    # Valores do vento a 1 UA, compartilhados pelo cálculo do empuxo.
    # São valores padrão até o serviço os atualizar com a amostra da missão.
    _v_sw_km_s = 400.0
    _n_earth_m3 = 7.3e6
    _t_e_earth_ev = 10.0
    _config_source = "default"
    _ecliptic_normal_i = np.array([0.0, 0.0, 1.0])

    @staticmethod
    def _tcc_orbital_basis(r_inercial):
        """Retorna a base O_o do TCC escrita em coordenadas inerciais.

        O TCC define z_o na direcao radial Sol-sonda. O eixo y_o e
        perpendicular simultaneamente a z_o e a normal da ecliptica; x_o
        completa a base dextrogira.
        """
        r_norm = np.linalg.norm(r_inercial)
        if not np.isfinite(r_norm) or r_norm <= 0.0:
            raise ValueError("A posicao da sonda deve ter norma positiva.")

        z_o = np.asarray(r_inercial, dtype=float) / r_norm
        y_raw = np.cross(ElectricSailDynamic._ecliptic_normal_i, z_o)
        y_norm = np.linalg.norm(y_raw)

        # Quando z_o e paralelo a z_i, a definicao geometrica do TCC nao
        # determina y_o: nao existe interseccao unica entre os planos.
        if not np.isfinite(y_norm) or y_norm <= 1e-12:
            raise ValueError(
                "A base orbital do TCC e indefinida quando r e paralelo "
                "a normal da ecliptica."
            )

        y_o = y_raw / y_norm
        x_o = np.cross(y_o, z_o)
        x_norm = np.linalg.norm(x_o)
        if not np.isfinite(x_norm) or x_norm <= 1e-12:
            raise ValueError("Nao foi possivel construir a base orbital do TCC.")
        x_o = x_o / x_norm

        # As colunas sao os eixos de O_o escritos em O_i.
        return np.column_stack((x_o, y_o, z_o))

    @classmethod
    def configure_from_csv(cls, start_time, csv_path=None):
        """Carrega a amostra mais próxima no mesmo dia e atualiza o vento a 1 UA.

        A tabela fornece temperatura em K, densidade em cm^-3 e velocidade em
        km/s; aqui são convertidas para eV, m^-3 e km/s, respectivamente.
        """
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
        dens_cm3 = float(row[4])  # densidade de prótons da tabela, usada como n no modelo
        vel_km_s = float(row[5])

        cls._v_sw_km_s = vel_km_s  # km/s
        cls._n_earth_m3 = dens_cm3 * 1e6  # cm^-3 -> m^-3
        cls._t_e_earth_ev = temp_k / 11604.5  # K -> eV
        cls._config_source = f"{year}-DOY{doy}-H{int(round(row[2]))}"

    @staticmethod
    def calculate_thrust_per_m(body):
        """
        Calcula a força de empuxo por metro de fio (σ_F) em N/m.

        Usa os valores de referência a 1 UA. A dependência radial aproximada
        do empuxo é aplicada em calculate_acceleration, conforme a eq. 73.
        """
        # O modelo interno usa SI: km/s -> m/s e eV -> J onde necessário.
        v_sw = ElectricSailDynamic._v_sw_km_s * 1e3
        n_earth = ElectricSailDynamic._n_earth_m3
        T_e_earth_eV = ElectricSailDynamic._t_e_earth_ev

        n = n_earth
        T_e_joules = T_e_earth_eV * E_CARGA_F

        # Parametros da vela
        V = body.V_volts  # tensão elétrica (V)
        r_w = body.r_w_km * 1000  # km -> m

        termo_raiz_ln = np.sqrt((EPSILON_0 * T_e_joules) / (n * E_CARGA_F**2))
        arg_ln = (2 / r_w) * termo_raiz_ln
        resultado_ln = np.log(arg_ln)
        expoente = (M_PROTON * v_sw**2 / (E_CARGA_F * V)) * resultado_ln
        denominador = E_CARGA_F * np.sqrt(np.exp(expoente) - 1)
        numerador = 6.18 * M_PROTON * v_sw**2 * np.sqrt(n * EPSILON_0 * T_e_joules)

        sigma_F = numerador / denominador

        return sigma_F

    @staticmethod
    def calculate_acceleration(body: ElectricSailProbe):
        """
        Calcula a aceleração da e-sail em coordenadas inerciais, em km/s^2.

        A massa propulsiva é calculada pelas equações 75–76; não é body.mass.
        """
        # Parâmetros da vela
        N = body.N
        L = body.L * 1000  # comprimento por fio: km -> m
        phi = body.phi
        theta = body.theta
        r_w = body.r_w_km * 1000  # raio do fio: km -> m
        V = body.V_volts  # tensão elétrica em volts (não é velocidade)

        # A integração usa km; as equações de empuxo usam metros.
        r_inercial = body.position * 1000  # km -> m
        
        # Distancia ao Sol
        r = np.linalg.norm(r_inercial)

        # Empuxo por metro calculado com valores de referência; a lei radial
        # (r_AU/r)^(7/6) é aplicada uma única vez abaixo.
        sigma_F_base = ElectricSailDynamic.calculate_thrust_per_m(body)

        magnitude_forca = (1/2) * N * L * sigma_F_base * (AU_M / r)**(7/6)
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)

        Fx_orbita = cos_phi * sin_theta * cos_theta
        Fy_orbita = -sin_phi * cos_phi * cos_theta**2
        Fz_orbita = cos_phi**2 * cos_theta**2 + 1
        
        F_vela_orbita = magnitude_forca * np.array([Fx_orbita, Fy_orbita, Fz_orbita])
        # Matriz de rotacao O_o -> O_i conforme a definicao geometrica do TCC
        matriz_rotacao = ElectricSailDynamic._tcc_orbital_basis(r_inercial)

        # rotacao
        F_vela_inercial = matriz_rotacao @ F_vela_orbita

        # A massa usada vem das eqs. 75-78, não de body.mass.
        # Massa do corpo da vela por metro do fio (eq 75)
        sigma_mb = 2 * K_T * BETA * N_EARTH * r_w * np.sqrt((2 * E_CARGA**3 * V**3) / M_ELETRON)
        
        # Massa dos fios por metro de fio (eq 76)
        sigma_mt = K_T * np.pi * RHO_W * r_w**2
        
        # Massa seca total por metro de fio
        sigma_m_vela = sigma_mb + sigma_mt
        
        L_total = N * L
        m_vela = sigma_m_vela * L_total
        
        # a = F/m
        a_vela_inercial = (F_vela_inercial / m_vela) * 0.001  # converte de m/s^2 para km/s^2

        return a_vela_inercial
