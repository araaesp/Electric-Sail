"""
Ferramentas auxiliares para definir a orientação (phi, theta) da vela elétrica
no referencial Local-Vertical/Local-Horizontal (LOF).

- r̂: direção radial do Sol para a sonda.
- t̂: direção tangencial (prograda) obtida pela projeção da velocidade no plano orbital.
- n̂: direção normal ao plano orbital (r̂ × t̂).

A força calculada em `ElectricSailDynamic` obedece à equação (102):

    F_LOF = prefator * [ T, N, R ] =
             prefator * [
                cos(phi) * sin(theta) * cos(theta),
                -sin(phi) * cos(phi) * cos^2(theta),
                cos^2(phi) * cos^2(theta) + 1
             ]

Este módulo oferece:

1) Construtores de base LOF e normalização de vetores.
2) Funções para converter (phi, theta) em direção de força.
3) Uma busca grosseira para encontrar (phi, theta) que melhor alinhem a força
   desejada (T, N, R) com um vetor alvo fornecido.
4) Heurísticas simples para fases típicas (espiralar para fora, para dentro,
   mudar inclinação, etc.).

Ainda não está conectado ao solver; serve como bloco isolado para análise
e integração futura.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable, Optional, Tuple

import numpy as np

AnglesDeg = Tuple[float, float]
Vector3 = np.ndarray


def _normalize(vec: Vector3) -> Vector3:
    vec = np.asarray(vec, dtype=float)
    norm = np.linalg.norm(vec)
    if norm == 0:
        raise ValueError("Não é possível normalizar vetor de norma zero.")
    return vec / norm


def compute_lof_basis(r_vec_km: Vector3, v_vec_km_s: Vector3) -> Tuple[Vector3, Vector3, Vector3]:
    """
    Retorna (t̂, n̂, r̂) em coordenadas inerciais dadas posição e velocidade em km.
    """
    r_hat = _normalize(r_vec_km)
    h_vec = np.cross(r_vec_km, v_vec_km_s)
    n_hat = _normalize(h_vec)
    t_hat = _normalize(np.cross(n_hat, r_hat))
    return t_hat, n_hat, r_hat


def sail_force_components(phi_rad: float, theta_rad: float) -> Vector3:
    """
    Componentes (T, N, R) da força na base LOF para ângulos em radianos.
    """
    cos_phi = np.cos(phi_rad)
    sin_phi = np.sin(phi_rad)
    cos_theta = np.cos(theta_rad)
    sin_theta = np.sin(theta_rad)

    t_comp = cos_phi * sin_theta * cos_theta
    n_comp = -sin_phi * cos_phi * cos_theta**2
    r_comp = cos_phi**2 * cos_theta**2 + 1.0
    return np.array([t_comp, n_comp, r_comp])


def sail_force_direction(phi_rad: float, theta_rad: float) -> Vector3:
    """
    Direção unitária da força (T, N, R) dada a orientação da vela.
    """
    comps = sail_force_components(phi_rad, theta_rad)
    return _normalize(comps)


def find_angles_for_direction(
    target_lof: Iterable[float],
    grid_step_deg: float = 0.5,
    theta_limits_deg: Tuple[float, float] = (0.5, 179.5),
    phi_limits_deg: Tuple[float, float] = (-89.5, 89.5),
) -> Tuple[AnglesDeg, float]:
    """
    Busca exaustiva (grade) que encontra (phi, theta) em graus cujo vetor
    de força no LOF se alinha ao máximo com `target_lof`.

    Retorna ((phi_deg, theta_deg), custo), onde custo = 1 - cos(erro).
    """
    target_dir = _normalize(np.asarray(target_lof, dtype=float))

    phi_values = np.deg2rad(
        np.arange(phi_limits_deg[0], phi_limits_deg[1] + grid_step_deg, grid_step_deg)
    )
    theta_values = np.deg2rad(
        np.arange(theta_limits_deg[0], theta_limits_deg[1] + grid_step_deg, grid_step_deg)
    )

    best_cost = np.inf
    best_angles = (0.0, 90.0)

    for phi in phi_values:
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)
        cos_phi_sq = cos_phi**2

        # Pré-calcula termos que dependem de phi para acelerar a varredura em theta.
        for theta in theta_values:
            cos_theta = np.cos(theta)
            sin_theta = np.sin(theta)

            t_comp = cos_phi * sin_theta * cos_theta
            n_comp = -sin_phi * cos_phi * cos_theta**2
            r_comp = cos_phi_sq * cos_theta**2 + 1.0

            vec = np.array([t_comp, n_comp, r_comp])
            vec_dir = vec / np.linalg.norm(vec)

            cost = 1.0 - float(np.dot(vec_dir, target_dir))
            if cost < best_cost:
                best_cost = cost
                best_angles = (np.rad2deg(phi), np.rad2deg(theta))

    return best_angles, best_cost


@dataclass
class AttitudeCommand:
    """
    Representa um comando de atitude desejado para a vela: direção alvo no LOF.
    """
    description: str
    target_lof: Vector3

    def solve_angles(self, **kwargs) -> Tuple[AnglesDeg, float]:
        """Resolve (phi, theta) para o alvo atual usando `find_angles_for_direction`."""
        return find_angles_for_direction(self.target_lof, **kwargs)


def heuristic_spiral_outward() -> AttitudeCommand:
    """
    Maximiza componente tangencial prograda mantendo força radial positiva.
    Aproximação: vetor alvo ~ [1, 0, 1].
    """
    target = _normalize(np.array([1.0, 0.0, 1.0]))
    return AttitudeCommand(
        description="Espiralar para fora (prograda)",
        target_lof=target,
    )


def heuristic_spiral_inward() -> AttitudeCommand:
    """
    Componente tangencial retrógrada para reduzir semieixo maior.
    """
    target = _normalize(np.array([-1.0, 0.0, 1.0]))
    return AttitudeCommand(
        description="Espiralar para dentro (retrógrada)",
        target_lof=target,
    )


def heuristic_plane_change(upwards: bool = True, tangential_weight: float = 0.1) -> AttitudeCommand:
    """
    Produz componente normal dominante para alterar inclinação.
    O parâmetro `tangential_weight` adiciona leve componente tangencial
    para evitar perda completa de controle de energia.
    """
    sign = 1.0 if upwards else -1.0
    target = _normalize(np.array([tangential_weight, sign, 1.0]))
    return AttitudeCommand(
        description=f"Alterar inclinação ({'subir' if upwards else 'descer'})",
        target_lof=target,
    )