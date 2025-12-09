"""
RFT cosmology parameter definitions.

This module defines a single dataclass RFTParameters that collects all
model and calibration parameters used in the RFT_UnifiedDeformation
background evolution and observable calculations.

All parameters have explicit units or are dimensionless. Numerical values
are chosen to produce a deformation strength Gamma_eff of order 10^-2 and
a BAO scale near 147 Mpc for the reference configuration, while keeping
the system numerically stable.
"""

from dataclasses import dataclass


@dataclass
class RFTParameters:
    # Fundamental constants (natural units for the core)
    c: float = 1.0          # speed of light
    hbar: float = 1.0       # reduced Planck constant
    G: float = 1.0          # Newton's constant in the core equations

    # Effective scalar masses (dimensionless in core units)
    m_Phi: float = 1.2
    m_Gamma: float = 0.9
    m_R: float = 0.8
    m_Xi: float = 0.7
    m_Psi: float = 0.6

    # Quadratic cross–couplings (dimensionless, symmetric)
    lambda_GammaR: float = 0.15
    lambda_PhiGamma: float = 0.12
    lambda_XiPsi: float = 0.10

    # Effective cosmological constant term (dimensionless in core units)
    Lambda_RFT: float = 1.0e-4

    # Metric couplings for conformal factor A(Phi,R) = 1 + a_Phi*Phi + a_R*R
    a_Phi: float = 0.05
    a_R: float = 0.05

    # Deformation of effective Planck mass M_*^2(Gamma) = M_Pl^2 (1 + xi_Gamma*Gamma)
    xi_Gamma: float = 0.20

    # Present–day (dimensionless) energy densities in core units
    # These play the role of Omega_b0 * (3 M_Pl^2 H0^2), etc., in the toy model.
    rho_b0: float = 0.30
    rho_gamma0: float = 5.0e-5

    # Initial conditions for the homogeneous fields at t = 0 (early time)
    Phi_0: float = 1.0e-3
    Gamma_0: float = 2.0e-3
    R_0: float = 1.0e-3
    Xi_0: float = 5.0e-4
    Psi_0: float = 2.0e-3

    # Initial time derivatives for the fields at t = 0
    Phi_dot_0: float = 0.0
    Gamma_dot_0: float = 0.0
    R_dot_0: float = 0.0
    Xi_dot_0: float = 0.0
    Psi_dot_0: float = 0.0

    # Initial scale factor and Hubble rate in core units at t = 0
    a_0: float = 1.0e-4
    H_0: float = 1.0e-2

    # Reference deformation scalar for calibration
    Gamma_ref: float = 1.9904e-2

    # BAO calibration: r_BAO(Gamma_ref) = r_BAO_ref_mpc
    r_BAO_ref_mpc: float = 147.0

    # Galaxy rotation calibration: a0(Gamma_ref) = a0_ref in (km/s)^2/kpc
    a0_ref: float = 2500.0

    # Collapse calibration: lambda_RFT(1 g, 1 µm; Gamma_ref) = lambda0_ref s^-1
    lambda0_ref: float = 1.0e-10

    # Reference mass and separation for collapse scaling
    m_ref_g: float = 1.0     # grams
    dx_ref_um: float = 1.0   # microns

    # Galaxy disk parameters for the example rotation curve
    M_disk_Msun: float = 1.2e11   # solar masses
    R_disk_kpc: float = 3.0       # kpc

    # Gravitational constant in astrophysical units for rotation curves
    # G_phys = 4.3009e-6 kpc (km/s)^2 / Msun
    G_kpc: float = 4.3009e-6      # kpc (km/s)^2 / Msun


# Default parameter set used in the examples and top–level functions.
DEFAULT_PARAMS = RFTParameters()
