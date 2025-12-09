"""
Minimal ΛCDM comparison utilities.

This module provides simple ΛCDM formulas for:
  - sound horizon at recombination,
  - acoustic scale ell_A,
  - Newtonian baryon–only rotation curves.

The goal is honest side–by–side comparison of characteristic scales,
not a full precision ΛCDM pipeline.
"""

from typing import Tuple

import numpy as np

from .parameters import RFTParameters, DEFAULT_PARAMS


def lcdm_rs_mpc(h: float = 0.68,
                Omega_b: float = 0.049,
                Omega_cdm: float = 0.262,
                Omega_gamma: float = 5.0e-5) -> float:
    """
    Approximate ΛCDM sound horizon at recombination in Mpc.

    Uses a standard fitting–formula level approximation for r_s
    at the drag epoch. This is not intended as a replacement for CAMB/CLASS,
    only as a reference scale.
    """
    # Eisenstein & Hu–style approximate formula for drag epoch sound horizon.
    # This is a widely used approximation.
    omega_b = Omega_b * h**2
    omega_m = (Omega_b + Omega_cdm) * h**2

    b1 = 0.313 * (omega_m**-0.419) * (1.0 + 0.607 * (omega_m**0.674))
    b2 = 0.238 * (omega_m**0.223)
    z_drag = 1291.0 * (omega_m**0.251) / (1.0 + 0.659 * (omega_m**0.828)) * (1.0 + b1 * (omega_b**b2))

    R_drag = 31.5 * omega_b * 1.0e3 / z_drag
    R_eq = 31.5 * omega_b * 1.0e3 / 3400.0
    k_eq = 0.0746 * omega_m

    rs = (2.0 / (3.0 * k_eq)) * np.sqrt(6.0 / R_eq) * np.log(
        (np.sqrt(1.0 + R_drag) + np.sqrt(R_drag + R_eq)) / (1.0 + np.sqrt(R_eq))
    )

    return float(rs)


def lcdm_ell_A(D_A_mpc: float, rs_mpc: float) -> float:
    """
    Acoustic scale ell_A in ΛCDM:
        ell_A = pi * D_A(z*) / r_s(z*)
    """
    if rs_mpc <= 0.0:
        return float("nan")
    return float(np.pi * D_A_mpc / rs_mpc)


def lcdm_rotation_curve_newtonian(r_kpc: np.ndarray,
                                  M_disk_Msun: float,
                                  R_disk_kpc: float,
                                  G_kpc: float) -> np.ndarray:
    """
    Pure Newtonian baryon–only rotation curve for an exponential disk.

    This uses the same baryon distribution as the RFT rotation curve for
    direct comparison. No dark matter halo is included here.
    """
    x = r_kpc / R_disk_kpc
    M_r = M_disk_Msun * (1.0 - (1.0 + x) * np.exp(-x))
    gN = G_kpc * M_r / (r_kpc**2)
    v_newton = np.sqrt(r_kpc * gN)
    return v_newton
