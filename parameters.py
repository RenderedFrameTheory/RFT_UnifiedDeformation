"""
RFT minimal parameter file for the unified deformation testbed.

Units: natural units with c = ħ = G = 1 for the dynamical core.
All alpha_i, beta_i are dimensionless effective couplings in render time τ.
Initial conditions are dimensionless field amplitudes at τ = 0,
corresponding to an early-universe state (e.g. t ~ 10^-35 s or η ~ 10^-33).

This file deliberately contains ONLY:
  - fundamental constants (set to 1 here),
  - alpha1..alpha5, beta1..beta5,
  - Phi_0, Gamma_0, R_0, Xi_0, Psi_0.
"""

# -------------------------------------------------------------------------
# Fundamental constants (set to 1 in this minimal model)
# -------------------------------------------------------------------------
c    = 1.0   # speed of light
hbar = 1.0   # reduced Planck constant
G    = 1.0   # Newton's constant

# -------------------------------------------------------------------------
# Effective RFT couplings for the 5-field coherence-flow system
# -------------------------------------------------------------------------
# Damping / self-couplings (dimensionless)
alpha1 = 0.80   # Φ damping rate
alpha2 = 0.60   # Γ damping rate
alpha3 = 0.50   # R damping rate
alpha4 = 0.40   # Ξ damping rate
alpha5 = 0.35   # Ψ damping rate

# Cross-couplings (dimensionless)
beta1  = 0.30   # Γ → Φ coupling
beta2  = 0.25   # Φ ↔ R coupling
beta3  = 0.20   # Γ ↔ Ξ coupling
beta4  = 0.18   # R ↔ Ψ coupling
beta5  = 0.15   # Ξ ↔ Ψ coupling

# -------------------------------------------------------------------------
# Initial conditions at early render time τ = 0
# -------------------------------------------------------------------------
Phi_0   = 1.0e-3   # initial render potential
Gamma_0 = 2.0e-3   # initial deformation field
R_0     = 1.0e-3   # initial geometric response
Xi_0    = 5.0e-4   # initial coherence field
Psi_0   = 2.0e-3   # initial observer-capacity field
