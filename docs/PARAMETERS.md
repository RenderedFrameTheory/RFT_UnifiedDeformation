# RFT Cosmology Parameters

This document lists the parameters defined in `rft/parameters.py`, their units, values and roles.

## 1. Fundamental constants

| Symbol | Name              | Code name | Value | Unit | Role |
|--------|-------------------|-----------|-------|------|------|
| c      | Speed of light    | c         | 1.0   | core units | Sets the light speed; absorbed into units. |
| ħ      | Reduced Planck    | hbar      | 1.0   | core units | Sets quantum scale; absorbed into units. |
| G      | Newton's constant | G         | 1.0   | core units | Used in background units; physical G_kpc used for rotation curves. |

## 2. Scalar masses

| Symbol  | Code name | Value | Unit        | Role |
|---------|-----------|-------|-------------|------|
| m_Φ     | m_Phi     | 1.2   | core mass   | Sets stiffness of Φ potential. |
| m_Γ     | m_Gamma   | 0.9   | core mass   | Controls deformation field response. |
| m_R     | m_R       | 0.8   | core mass   | Controls curvature–response field. |
| m_Ξ     | m_Xi      | 0.7   | core mass   | Controls coherence field dynamics. |
| m_Ψ     | m_Psi     | 0.6   | core mass   | Controls observer–capacity field. |

## 3. Cross–couplings

| Coupling             | Code name       | Value | Unit       | Role |
|----------------------|-----------------|-------|------------|------|
| λ_{ΓR}               | lambda_GammaR   | 0.15  | core units | Couples Γ to R; mixes deformation and curvature response. |
| λ_{ΦΓ}               | lambda_PhiGamma | 0.12  | core units | Couples Φ to Γ; links render potential to deformation. |
| λ_{ΞΨ}               | lambda_XiPsi    | 0.10  | core units | Couples Ξ to Ψ; links coherence and observer fields. |

## 4. Constant term

| Symbol        | Code name   | Value   | Unit       | Role |
|---------------|-------------|---------|------------|------|
| Λ_{RFT}       | Lambda_RFT  | 1.0e-4  | core units | Effective cosmological constant in the scalar sector. |

## 5. Metric couplings

| Symbol   | Code name | Value | Unit       | Role |
|----------|-----------|-------|------------|------|
| a_Φ      | a_Phi     | 0.05  | core units | Controls how Φ deforms the physical metric seen by matter. |
| a_R      | a_R       | 0.05  | core units | Controls how R deforms the physical metric seen by matter. |
| ξ_Γ      | xi_Gamma  | 0.20  | core units | Controls how Γ modifies the effective Planck mass. |

## 6. Present–day densities (Einstein frame, core units)

| Quantity       | Code name  | Value   | Role |
|----------------|------------|---------|------|
| ρ_{b0}         | rho_b0     | 0.30    | Effective baryon density scale. |
| ρ_{γ0}         | rho_gamma0 | 5.0e-5  | Effective radiation density scale. |

These values are chosen to produce a matter–radiation mix qualitatively similar to standard cosmology in the toy units of the model.

## 7. Initial conditions

At t = 0:

| Quantity     | Code name   | Value   | Role |
|--------------|-------------|---------|------|
| a(0)         | a_0         | 1.0e-4  | Initial scale factor. |
| H(0)         | H_0         | 1.0e-2  | Initial Hubble rate. |
| Φ(0)         | Phi_0       | 1.0e-3  | Initial render potential. |
| Γ(0)         | Gamma_0     | 2.0e-3  | Initial deformation field. |
| R(0)         | R_0         | 1.0e-3  | Initial curvature–response field. |
| Ξ(0)         | Xi_0        | 5.0e-4  | Initial coherence field. |
| Ψ(0)         | Psi_0       | 2.0e-3  | Initial observer–capacity field. |
| Φ̇(0)        | Phi_dot_0   | 0.0     | Initial Φ velocity. |
| Γ̇(0)        | Gamma_dot_0 | 0.0     | Initial Γ velocity. |
| Ṙ(0)        | R_dot_0     | 0.0     | Initial R velocity. |
| Ξ̇(0)        | Xi_dot_0    | 0.0     | Initial Ξ velocity. |
| Ψ̇(0)        | Psi_dot_0   | 0.0     | Initial Ψ velocity. |

These values define a stable starting point for the background integration without fine tuning.

## 8. Deformation and calibration parameters

| Symbol                 | Code name        | Value      | Role |
|------------------------|------------------|------------|------|
| Γ_{ref}                | Gamma_ref        | 1.9904e-2  | Reference deformation scalar for calibration. |
| r_{BAO,ref}            | r_BAO_ref_mpc    | 147.0      | BAO scale at Γ_ref in Mpc. |
| a_{0,ref}              | a0_ref           | 2500.0     | Galaxy acceleration scale at Γ_ref in (km/s)^2/kpc. |
| λ_{0,ref}              | lambda0_ref      | 1.0e-10    | Collapse rate at Γ_ref, m_ref and dx_ref in s^-1. |
| m_ref                  | m_ref_g          | 1.0        | Reference mass for collapse, grams. |
| Δx_ref                 | dx_ref_um        | 1.0        | Reference separation for collapse, microns. |

## 9. Galaxy disk and gravitational constant (astrophysical units)

| Quantity          | Code name      | Value   | Unit                      |
|-------------------|----------------|---------|---------------------------|
| M_disk            | M_disk_Msun    | 1.2e11  | Solar masses              |
| R_disk            | R_disk_kpc     | 3.0     | kpc                       |
| G_phys            | G_kpc          | 4.3009e-6 | kpc (km/s)^2 / Msun     |

These parameters define the Milky Way–like baryonic disk used in the example rotation curves.
