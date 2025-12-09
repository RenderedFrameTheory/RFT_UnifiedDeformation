"""
Numerical evolution of the RFT background system.

This module provides functions to integrate the FRW+five–scalar system,
save the background solution to a NPZ file, and reload it for use in
observable calculations.
"""

from typing import Tuple

import numpy as np
from scipy.integrate import solve_ivp

from .parameters import RFTParameters, DEFAULT_PARAMS
from .cosmology_model import (
    rft_background_rhs,
    initial_background_state,
    compute_energy_densities,
    compute_gamma_eff,
    A_conformal,
)


def run_background_evolution(
    t_start: float = 0.0,
    t_end: float = 1.0e4,
    n_steps: int = 5000,
    params: RFTParameters = DEFAULT_PARAMS,
) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Integrate the RFT background system from t_start to t_end.

    Returns
    -------
    t_arr : ndarray
        Time samples.
    y_arr : ndarray
        Background state samples with shape (n_state, n_time).
    Gamma_eff : float
        Effective deformation scalar computed from the Gamma(t) trajectory.
    """
    y0 = initial_background_state(params)

    t_eval = np.linspace(t_start, t_end, n_steps)

    def rhs(t, y):
        return rft_background_rhs(t, y, params)

    sol = solve_ivp(
        rhs,
        t_span=(t_start, t_end),
        y0=y0,
        t_eval=t_eval,
        method="RK45",
        rtol=1.0e-7,
        atol=1.0e-9,
    )

    if not sol.success:
        raise RuntimeError(f"RFT background evolution failed: {sol.message}")

    t_arr = sol.t
    y_arr = sol.y  # shape (12, n_steps)

    # Extract a(t), Gamma(t) and rho_b(t) for Gamma_eff calculation
    a_arr = y_arr[0]
    Gamma_arr = y_arr[4]

    rho_b_arr = []
    for k in range(t_arr.size):
        rho_b, rho_gamma, rho_rft, rho_tot = compute_energy_densities(y_arr[:, k], params)
        rho_b_arr.append(rho_b)
    rho_b_arr = np.array(rho_b_arr, dtype=float)

    Gamma_eff = compute_gamma_eff(t_arr, Gamma_arr, a_arr, rho_b_arr)

    return t_arr, y_arr, Gamma_eff


def save_background_solution(path: str,
                             t_arr: np.ndarray,
                             y_arr: np.ndarray,
                             Gamma_eff: float,
                             params: RFTParameters = DEFAULT_PARAMS) -> None:
    """
    Save the background solution and Gamma_eff to a NPZ file.

    The file contains:
        t           : time samples
        y           : state array (12 x N)
        Gamma_eff   : scalar deformation parameter
        a           : scale factor samples
        Phi, Gamma, R, Xi, Psi
        Phi_dot, Gamma_dot, R_dot, Xi_dot, Psi_dot
        rho_b, rho_gamma, rho_rft, rho_tot
    """
    a_arr = y_arr[0]
    Phi_arr = y_arr[2]
    Gamma_arr = y_arr[4]
    R_arr = y_arr[6]
    Xi_arr = y_arr[8]
    Psi_arr = y_arr[10]

    Phi_dot_arr = y_arr[3]
    Gamma_dot_arr = y_arr[5]
    R_dot_arr = y_arr[7]
    Xi_dot_arr = y_arr[9]
    Psi_dot_arr = y_arr[11]

    rho_b_arr = []
    rho_gamma_arr = []
    rho_rft_arr = []
    rho_tot_arr = []

    for k in range(t_arr.size):
        rho_b, rho_gamma, rho_rft, rho_tot = compute_energy_densities(y_arr[:, k], params)
        rho_b_arr.append(rho_b)
        rho_gamma_arr.append(rho_gamma)
        rho_rft_arr.append(rho_rft)
        rho_tot_arr.append(rho_tot)

    rho_b_arr = np.array(rho_b_arr, dtype=float)
    rho_gamma_arr = np.array(rho_gamma_arr, dtype=float)
    rho_rft_arr = np.array(rho_rft_arr, dtype=float)
    rho_tot_arr = np.array(rho_tot_arr, dtype=float)

    np.savez(
        path,
        t=t_arr,
        y=y_arr,
        Gamma_eff=Gamma_eff,
        a=a_arr,
        Phi=Phi_arr,
        Gamma=Gamma_arr,
        R=R_arr,
        Xi=Xi_arr,
        Psi=Psi_arr,
        Phi_dot=Phi_dot_arr,
        Gamma_dot=Gamma_dot_arr,
        R_dot=R_dot_arr,
        Xi_dot=Xi_dot_arr,
        Psi_dot=Psi_dot_arr,
        rho_b=rho_b_arr,
        rho_gamma=rho_gamma_arr,
        rho_rft=rho_rft_arr,
        rho_tot=rho_tot_arr,
    )


def load_background_solution(path: str) -> Tuple[np.ndarray, np.ndarray, float, dict]:
    """
    Load a saved background solution from a NPZ file.

    Returns
    -------
    t_arr : ndarray
    y_arr : ndarray
    Gamma_eff : float
    extras : dict
        Dictionary with additional arrays (a, Phi, Gamma, etc.).
    """
    data = np.load(path)
    t_arr = data["t"]
    y_arr = data["y"]
    Gamma_eff = float(np.asarray(data["Gamma_eff"]).squeeze())

    extras_keys = [
        "a", "Phi", "Gamma", "R", "Xi", "Psi",
        "Phi_dot", "Gamma_dot", "R_dot", "Xi_dot", "Psi_dot",
        "rho_b", "rho_gamma", "rho_rft", "rho_tot",
    ]
    extras = {key: data[key] for key in extras_keys if key in data.files}
    return t_arr, y_arr, Gamma_eff, extras
