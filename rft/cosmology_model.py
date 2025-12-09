"""
Core RFT cosmology model: action ingredients and FRW background equations.

This module implements the background–level equations for the five scalar
fields Φ, Γ, R, Ξ, Ψ and the scale factor a(t). It also provides the
definition and numerical evaluation of the deformation scalar Gamma_eff.

The implementation follows a scalar–tensor cosmology with:
  - physical metric  g_{μν}  in FRW form,
  - matter coupled to a conformally related metric
        tilde{g}_{μν} = A^2(Φ,R) g_{μν},
  - effective Planck mass M_*^2(Γ) = M_Pl^2 (1 + xi_Gamma * Γ).

The background state vector is
    y = [a, H,
         Φ, Φdot,
         Γ, Γdot,
         R, Rdot,
         Ξ, Ξdot,
         Ψ, Ψdot].
"""

from typing import Tuple, Callable

import numpy as np

from .parameters import RFTParameters


def A_conformal(Phi: float, R_field: float, p: RFTParameters) -> float:
    """
    Conformal factor A(Φ,R) that links the Einstein–frame metric g_{μν}
    to the physical metric seen by matter:
        tilde{g}_{μν} = A^2 g_{μν}
    """
    return 1.0 + p.a_Phi * Phi + p.a_R * R_field


def M_star_sq(Gamma: float, p: RFTParameters) -> float:
    """
    Effective Planck mass squared:
        M_*^2(Γ) = M_Pl^2 (1 + xi_Gamma * Γ)
    In this minimal implementation the overall M_Pl^2 factor is absorbed
    into the choice of units, so we keep only (1 + xi_Gamma*Gamma).
    """
    return 1.0 + p.xi_Gamma * Gamma


def potential(Phi: float, Gamma: float, R_field: float,
              Xi: float, Psi: float, p: RFTParameters) -> float:
    """
    Scalar potential:
        V = 1/2 m_Phi^2 Φ^2 + 1/2 m_Gamma^2 Γ^2 + 1/2 m_R^2 R^2
          + 1/2 m_Xi^2 Ξ^2 + 1/2 m_Psi^2 Ψ^2
          + lambda_GammaR Γ R + lambda_PhiGamma Φ Γ + lambda_XiPsi Ξ Ψ
          + Lambda_RFT
    """
    V_quad = (
        0.5 * p.m_Phi**2 * Phi**2 +
        0.5 * p.m_Gamma**2 * Gamma**2 +
        0.5 * p.m_R**2 * R_field**2 +
        0.5 * p.m_Xi**2 * Xi**2 +
        0.5 * p.m_Psi**2 * Psi**2
    )
    V_mix = (
        p.lambda_GammaR * Gamma * R_field +
        p.lambda_PhiGamma * Phi * Gamma +
        p.lambda_XiPsi * Xi * Psi
    )
    return V_quad + V_mix + p.Lambda_RFT


def dV_dPhi(Phi: float, Gamma: float, R_field: float,
            Xi: float, Psi: float, p: RFTParameters) -> float:
    return p.m_Phi**2 * Phi + p.lambda_PhiGamma * Gamma


def dV_dGamma(Phi: float, Gamma: float, R_field: float,
              Xi: float, Psi: float, p: RFTParameters) -> float:
    return p.m_Gamma**2 * Gamma + p.lambda_GammaR * R_field + p.lambda_PhiGamma * Phi


def dV_dR(Phi: float, Gamma: float, R_field: float,
          Xi: float, Psi: float, p: RFTParameters) -> float:
    return p.m_R**2 * R_field + p.lambda_GammaR * Gamma


def dV_dXi(Phi: float, Gamma: float, R_field: float,
           Xi: float, Psi: float, p: RFTParameters) -> float:
    return p.m_Xi**2 * Xi + p.lambda_XiPsi * Psi


def dV_dPsi(Phi: float, Gamma: float, R_field: float,
            Xi: float, Psi: float, p: RFTParameters) -> float:
    return p.m_Psi**2 * Psi + p.lambda_XiPsi * Xi


def energy_densities(a: float,
                     Phi: float, Phi_dot: float,
                     Gamma: float, Gamma_dot: float,
                     R_field: float, R_dot: float,
                     Xi: float, Xi_dot: float,
                     Psi: float, Psi_dot: float,
                     p: RFTParameters) -> Tuple[float, float, float, float]:
    """
    Compute Einstein–frame energy densities for baryons, radiation and
    the RFT scalar sector.

    The conformal factor A(Φ,R) multiplies physical densities as
        rho_E = A^4 rho_phys.
    """
    A = A_conformal(Phi, R_field, p)

    rho_b = A**4 * p.rho_b0 / (a**3)
    rho_gamma = A**4 * p.rho_gamma0 / (a**4)

    kinetic = 0.5 * (
        Phi_dot**2 + Gamma_dot**2 + R_dot**2 + Xi_dot**2 + Psi_dot**2
    )
    V = potential(Phi, Gamma, R_field, Xi, Psi, p)
    rho_rft = kinetic + V

    rho_tot = rho_b + rho_gamma + rho_rft
    return rho_b, rho_gamma, rho_rft, rho_tot


def rft_background_rhs(t: float,
                       y: np.ndarray,
                       p: RFTParameters) -> np.ndarray:
    """
    Right–hand side of the RFT background ODE system.

    State vector:
        y = [a, H,
             Φ, Φdot,
             Γ, Γdot,
             R, Rdot,
             Ξ, Ξdot,
             Ψ, Ψdot]
    """
    a, H, Phi, Phi_dot, Gamma, Gamma_dot, R_field, R_dot, Xi, Xi_dot, Psi, Psi_dot = y

    # Energy densities
    rho_b, rho_gamma, rho_rft, rho_tot = energy_densities(
        a, Phi, Phi_dot, Gamma, Gamma_dot, R_field, R_dot, Xi, Xi_dot, Psi, Psi_dot, p
    )

    # Effective Planck mass and its derivative
    Mstar2 = M_star_sq(Gamma, p)
    # In this simplified background we treat dot(M_*^2) as xi_Gamma * Gamma_dot
    Mstar2_dot = p.xi_Gamma * Gamma_dot

    # Ricci scalar for flat FRW: R = 6(2H^2 + dot{H})
    # We will use an approximate expression for dot{H} from the acceleration equation.
    # First compute kinetic term:
    kinetic = 0.5 * (
        Phi_dot**2 + Gamma_dot**2 + R_dot**2 + Xi_dot**2 + Psi_dot**2
    )

    # Approximate Raychaudhuri equation (Einstein frame, simplified)
    # -2 M_*^2 dot{H} ≈ rho_b + 4/3 rho_gamma + 2 kinetic
    # ignoring higher order terms in M_* variation.
    numerator = rho_b + (4.0 / 3.0) * rho_gamma + 2.0 * kinetic
    H_dot = -0.5 * numerator / max(Mstar2, 1.0e-10)

    # Ricci scalar in this approximation
    R_ricci = 6.0 * (2.0 * H**2 + H_dot)

    # Scalar equations of motion
    # Source from matter trace for Φ and R (baryons dominate)
    T_m = -rho_b  # pressureless baryons
    S_Phi = -p.a_Phi * T_m
    S_R = -p.a_R * T_m

    Phi_ddot = -3.0 * H * Phi_dot - dV_dPhi(Phi, Gamma, R_field, Xi, Psi, p) + S_Phi
    Gamma_ddot = (
        -3.0 * H * Gamma_dot
        - dV_dGamma(Phi, Gamma, R_field, Xi, Psi, p)
        + 0.5 * p.xi_Gamma * R_ricci
    )
    R_ddot = -3.0 * H * R_dot - dV_dR(Phi, Gamma, R_field, Xi, Psi, p) + S_R
    Xi_ddot = -3.0 * H * Xi_dot - dV_dXi(Phi, Gamma, R_field, Xi, Psi, p)
    Psi_ddot = -3.0 * H * Psi_dot - dV_dPsi(Phi, Gamma, R_field, Xi, Psi, p)

    # Scale factor evolution: dot{a} = a H
    a_dot = a * H

    # Assemble derivative vector
    return np.array([
        a_dot,
        H_dot,
        Phi_dot,
        Phi_ddot,
        Gamma_dot,
        Gamma_ddot,
        R_dot,
        R_ddot,
        Xi_dot,
        Xi_ddot,
        Psi_dot,
        Psi_ddot,
    ], dtype=float)


def initial_background_state(p: RFTParameters) -> np.ndarray:
    """
    Construct the initial state vector y(t=0) from the parameter set.
    """
    y0 = np.array([
        p.a_0,
        p.H_0,
        p.Phi_0,
        p.Phi_dot_0,
        p.Gamma_0,
        p.Gamma_dot_0,
        p.R_0,
        p.R_dot_0,
        p.Xi_0,
        p.Xi_dot_0,
        p.Psi_0,
        p.Psi_dot_0,
    ], dtype=float)
    return y0


def compute_energy_densities(y: np.ndarray,
                             p: RFTParameters) -> Tuple[float, float, float, float]:
    """
    Convenience wrapper to get Einstein–frame energy densities from a state y.
    """
    a, H, Phi, Phi_dot, Gamma, Gamma_dot, R_field, R_dot, Xi, Xi_dot, Psi, Psi_dot = y
    return energy_densities(
        a, Phi, Phi_dot, Gamma, Gamma_dot, R_field, R_dot, Xi, Xi_dot, Psi, Psi_dot, p
    )


def compute_gamma_eff(t_arr: np.ndarray,
                      Gamma_arr: np.ndarray,
                      a_arr: np.ndarray,
                      rho_b_arr: np.ndarray) -> float:
    """
    Compute the effective deformation scalar Gamma_eff as a weighted average
    of Gamma(t) over the background evolution.

    Definition:
        Gamma_eff = [∫ Gamma(η) W(η) dη] / [∫ W(η) dη]
    with conformal time η and weight W ∝ a^2 rho_b.

    In terms of cosmic time t we approximate:
        dη = dt / a,
        W(t) ∝ a^2(t) rho_b(t),
    so the numerator is ∫ Gamma(t) a(t) rho_b(t) dt and the denominator
    is ∫ a(t) rho_b(t) dt.

    t_arr is assumed to be strictly increasing.
    """
    t = np.asarray(t_arr)
    Gamma = np.asarray(Gamma_arr)
    a = np.asarray(a_arr)
    rho_b = np.asarray(rho_b_arr)

    W = a * rho_b

    dt = np.diff(t)
    # Trapezoidal integration
    W_mid = 0.5 * (W[1:] + W[:-1])
    Gamma_mid = 0.5 * (Gamma[1:] + Gamma[:-1])

    num = np.sum(Gamma_mid * W_mid * dt)
    den = np.sum(W_mid * dt)

    if den <= 0.0:
        return float("nan")

    return float(num / den)
