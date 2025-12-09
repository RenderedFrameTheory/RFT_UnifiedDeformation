import numpy as np
from scipy.integrate import solve_ivp

from parameters import (
    alpha1, alpha2, alpha3, alpha4, alpha5,
    beta1, beta2, beta3, beta4, beta5,
    Phi_0, Gamma_0, R_0, Xi_0, Psi_0,
)

# Render-time interval (dimensionless)
TAU_START = 0.0
TAU_END   = 1.0


def rft_rhs(tau, y):
    """
    RFT 5-field coherence-flow system in render time τ.

    State vector:
      y = [Phi, Gamma, R, Xi, Psi]

    Equations:
      dPhi/dτ   = -α1 Φ + β1 Γ - β2 R
      dGamma/dτ = -α2 Γ + β2 Φ - β3 Xi
      dR/dτ     = -α3 R + β3 Γ - β4 Psi
      dXi/dτ    = -α4 Xi + β4 R - β5 Phi
      dPsi/dτ   = -α5 Psi + β5 Xi
    """
    Phi, Gamma, R, Xi, Psi = y

    dPhi   = -alpha1 * Phi   + beta1 * Gamma - beta2 * R
    dGamma = -alpha2 * Gamma + beta2 * Phi   - beta3 * Xi
    dR     = -alpha3 * R     + beta3 * Gamma - beta4 * Psi
    dXi    = -alpha4 * Xi    + beta4 * R     - beta5 * Phi
    dPsi   = -alpha5 * Psi   + beta5 * Xi

    return np.array([dPhi, dGamma, dR, dXi, dPsi], dtype=float)


def compute_Gamma_eff(t_grid, Gamma_arr, weight_fn=None):
    """
    Compute Γ_eff as a weighted time average over τ.

    Γ_eff = [∫ W(τ) Γ(τ) dτ] / [∫ W(τ) dτ]
    """
    t = np.asarray(t_grid)
    G = np.asarray(Gamma_arr)

    if weight_fn is None:
        W = np.ones_like(t)
    else:
        W = weight_fn(t)

    dt = np.diff(t)
    G_mid = 0.5 * (G[1:] * W[1:] + G[:-1] * W[:-1])

    num = np.sum(G_mid * dt)
    den = np.sum(0.5 * (W[1:] + W[:-1]) * dt)

    return num / den


def run_evolution(
    n_steps: int = 2000,
    save_file: str = "rft_evolution.npz",
):
    """
    Integrate the RFT coherence system from TAU_START to TAU_END
    and save the result to an .npz file for plotting.
    """
    y0 = np.array([Phi_0, Gamma_0, R_0, Xi_0, Psi_0], dtype=float)

    t_eval = np.linspace(TAU_START, TAU_END, n_steps)

    sol = solve_ivp(
        fun=rft_rhs,
        t_span=(TAU_START, TAU_END),
        y0=y0,
        t_eval=t_eval,
        method="RK45",
        rtol=1e-8,
        atol=1e-10,
    )

    if not sol.success:
        raise RuntimeError(f"RFT evolution failed: {sol.message}")

    Phi   = sol.y[0]
    Gamma = sol.y[1]
    R     = sol.y[2]
    Xi    = sol.y[3]
    Psi   = sol.y[4]

    Gamma_eff = compute_Gamma_eff(sol.t, Gamma)

    np.savez(
        save_file,
        tau=sol.t,
        Phi=Phi,
        Gamma=Gamma,
        R=R,
        Xi=Xi,
        Psi=Psi,
        Gamma_eff=Gamma_eff,
    )

    print(f"[evolution] Saved evolution to {save_file}")
    print(f"[evolution] Γ_eff (time-avg) = {Gamma_eff:.6e}")

    return sol, Gamma_eff


if __name__ == "__main__":
    run_evolution()
