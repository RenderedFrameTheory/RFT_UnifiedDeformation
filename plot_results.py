import argparse
import math
import pathlib

import matplotlib.pyplot as plt
import numpy as np

# ------------------------------------------------------------
# GLOBAL REFERENCE VALUES (cosmology + RFT calibration)
# ------------------------------------------------------------

# Reference deformation strength for calibration
GAMMA_REF = 1.9904e-2

# Target BAO scale at reference Gamma (in Mpc)
R_BAO_REF_MPC = 147.0

# Compute BAO scale factor S_BAO so that r_BAO(GAMMA_REF) = 147 Mpc
R_BAO_DIMLESS_REF = math.pi / math.sqrt(GAMMA_REF)
S_BAO_MPC = R_BAO_REF_MPC / R_BAO_DIMLESS_REF

# Reference CMB peak positions at reference Gamma (Planck-like)
ELL_REF = np.array([220.0, 540.0, 800.0, 1100.0])

# ------------------------------------------------------------
# GALAXY ROTATION PARAMETERS (baryons only, RFT-modified gravity)
# ------------------------------------------------------------

# Gravitational constant in astrophysical units:
# G = 4.3009e-6 kpc (km/s)^2 / Msun
G_KPC = 4.3009e-6

# Exponential disk parameters (Milky Way-like, baryons only)
M_DISK_MSUN = 1.2e11   # total baryonic mass in solar masses
R_DISK_KPC  = 3.0      # disk scale length in kpc

# RFT acceleration scale a0(Gamma_eff) = A0_STAR * Gamma_eff
# A0_STAR chosen so that a0(GAMMA_REF) = A0_REF gives ~220-240 km/s flat tail
A0_REF = 2500.0                          # (km/s)^2 / kpc (effective)
A0_STAR = A0_REF / GAMMA_REF             # scale factor

# ------------------------------------------------------------
# COLLAPSE / COHERENCE PARAMETERS
# ------------------------------------------------------------

# Reference mass and separation
M_REF_G   = 1.0       # grams
DX_REF_UM = 1.0       # microns

# Reference collapse rate at Gamma_ref for m=1 g, dx=1 μm
LAMBDA0_REF = 1.0e-10   # s^-1

# ------------------------------------------------------------
# LOAD EVOLUTION
# ------------------------------------------------------------

def load_gamma_eff(path: str = "rft_evolution.npz") -> float:
    """
    Load and validate the scalar ``Gamma_eff`` value from an evolution archive.

    Parameters
    ----------
    path : str
        Path to the ``.npz`` file containing ``Gamma_eff``.

    Returns
    -------
    float
        Scalar ``Gamma_eff`` value.
    """
    file_path = pathlib.Path(path)
    if not file_path.exists():
        raise FileNotFoundError(
            f"Evolution file not found at '{file_path}'. "
            "Run the simulation pipeline first or specify --input."
        )

    with np.load(file_path) as data:
        if "Gamma_eff" not in data:
            raise KeyError(
                "Expected key 'Gamma_eff' in evolution file; available keys: "
                f"{list(data.keys())}"
            )

        gamma_raw = np.asarray(data["Gamma_eff"]).squeeze()

    try:
        gamma_val = float(gamma_raw)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "Gamma_eff must be convertible to a scalar float; "
            f"received type {type(gamma_raw)}"
        ) from exc

    if not math.isfinite(gamma_val) or gamma_val <= 0:
        raise ValueError(
            "Gamma_eff must be a finite, positive scalar for plotting; "
            f"received {gamma_val}"
        )

    return gamma_val


# ------------------------------------------------------------
# BAO SCALE & CMB PEAKS (Γ_eff-driven)
# ------------------------------------------------------------

def bao_scale_mpc(Gamma_eff: float) -> float:
    """
    Physical BAO scale r_BAO(Gamma_eff) in Mpc.

      r_dimless = π / sqrt(Gamma_eff)
      r_phys    = S_BAO_MPC * r_dimless

    S_BAO_MPC is fixed such that r_BAO(GAMMA_REF) ≈ 147 Mpc.
    """
    if Gamma_eff <= 0:
        return np.nan
    r_dimless = math.pi / math.sqrt(Gamma_eff)
    return S_BAO_MPC * r_dimless


def cmb_peak_locations(Gamma_eff: float) -> np.ndarray:
    """
    Scale reference peak locations with Gamma_eff:

      ell_n(Gamma) = ell_n_ref * sqrt(Gamma / GAMMA_REF)

    For Gamma_eff = GAMMA_REF, we recover ell_n_ref.
    """
    if Gamma_eff <= 0:
        return np.nan * np.ones_like(ELL_REF)
    scale = math.sqrt(Gamma_eff / GAMMA_REF)
    return ELL_REF * scale


def cmb_tt_template(Gamma_eff: float, ax=None):
    """
    Simple CMB TT template:
      - Gaussian peaks around ell_n(Gamma_eff)
      - Peak amplitudes/widths mildly modulated by Gamma_eff.
    """
    if ax is None:
        fig, ax = plt.subplots()

    ell = np.arange(2, 2501)
    ell_peaks = cmb_peak_locations(Gamma_eff)

    # Base amplitudes and widths
    base_amps = np.array([1.0, 0.65, 0.4, 0.25])
    base_width = 50.0

    # Modulate with Gamma_eff
    amp_scale   = 1.0 + 0.5 * Gamma_eff
    width_scale = 1.0 + 2.0 * Gamma_eff

    amps = base_amps * amp_scale
    width = base_width * width_scale

    Cl = np.zeros_like(ell, dtype=float)
    for A, ell0 in zip(amps, ell_peaks):
        Cl += A * np.exp(-0.5 * ((ell - ell0) / width) ** 2)

    D_ell = ell * (ell + 1) * Cl / (2.0 * math.pi)

    ax.plot(ell, D_ell, label="RFT TT template")
    for ell0 in ell_peaks:
        ax.axvline(ell0, linestyle="--", alpha=0.3)

    ax.set_xlabel(r"Multipole $\ell$")
    ax.set_ylabel(r"$\ell(\ell+1)C_\ell / 2\pi$ (arb. units)")
    ax.set_title("RFT CMB TT Power Spectrum (Γ_eff–tied peaks)")
    ax.set_xlim(2, 2500)
    ax.grid(True)
    ax.legend()

    return ax


def bao_plot(Gamma_eff: float, ax=None):
    """
    Plot matter power spectrum P(k) with a BAO bump tied to r_BAO(Gamma_eff).

    Template:
      P(k) ~ k^n_s T(k)^2 [1 + BAO_wiggle(k; k_BAO(Gamma_eff))].
    """
    if ax is None:
        fig, ax = plt.subplots()

    r_bao = bao_scale_mpc(Gamma_eff)   # in Mpc
    k_bao = 2.0 * math.pi / r_bao      # ~ h/Mpc (set h≈1 in this crude template)

    k = np.logspace(-4, 0, 400)        # h/Mpc

    n_s = 0.96
    P_prim = k ** n_s

    # Simple transfer function
    k_eq = 0.015
    T_k = 1.0 / (1.0 + (k / k_eq) ** 2)

    # BAO wiggle centred on k_bao
    A_bao = 0.1
    sigma_bao = 0.5 * k_bao
    bao_envelope = np.exp(-0.5 * ((k - k_bao) / sigma_bao) ** 2)
    bao_mod = 1.0 + A_bao * bao_envelope * np.sin((k / k_bao) * math.pi)

    P_k = P_prim * T_k**2 * bao_mod

    ax.loglog(k, P_k, label="RFT P(k) template")
    ax.axvline(k_bao, linestyle="--", alpha=0.5,
               label=fr"BAO scale ~ {r_bao:.1f} Mpc")

    ax.set_xlabel(r"$k\ (h/\mathrm{Mpc})$")
    ax.set_ylabel(r"$P(k)$ (arb. units)")
    ax.set_title("RFT Matter Power Spectrum with BAO (Γ_eff–tied)")
    ax.legend()
    ax.grid(True, which="both", ls=":")

    return ax


# ------------------------------------------------------------
# GALAXY ROTATION CURVE (RFT gravity, baryons only)
# ------------------------------------------------------------

def disk_mass_enclosed(r_kpc: np.ndarray,
                       M_disk: float = M_DISK_MSUN,
                       R_d: float = R_DISK_KPC) -> np.ndarray:
    """
    Exponential disk mass profile:
      M(<r) = M_disk * [1 - (1 + r/R_d) * exp(-r/R_d)]
    """
    x = r_kpc / R_d
    return M_disk * (1.0 - (1.0 + x) * np.exp(-x))


def rotation_curve_rft(Gamma_eff: float, ax=None):
    """
    Compute and plot RFT rotation curve for a Milky Way–like disk:

      - Baryonic exponential disk only.
      - Gravity modified via MOND-like 'simple μ' law:
          g_RFT = 0.5 * (g_N + sqrt(g_N^2 + 4*g_N*a0(Gamma_eff))).

      - a0(Gamma_eff) = A0_STAR * Gamma_eff, with A0_STAR calibrated
        so that for Gamma_ref, v_rot ~ 220–240 km/s from 10–40 kpc.
    """
    if ax is None:
        fig, ax = plt.subplots()

    # Radius grid (kpc)
    r = np.linspace(0.1, 100.0, 400)

    # Baryonic Newtonian acceleration
    M_r = disk_mass_enclosed(r)
    gN = G_KPC * M_r / r**2   # (km/s)^2 / kpc

    # RFT acceleration scale
    a0_eff = A0_STAR * Gamma_eff   # same units as gN

    # MOND-like solution with 'simple μ':
    #   g μ(g/a0) = gN, μ(x) = x/(1+x)
    # → g = 0.5 * [gN + sqrt(gN^2 + 4 gN a0)]
    gRFT = 0.5 * (gN + np.sqrt(gN**2 + 4.0 * gN * a0_eff))

    # Circular speeds
    v_RFT = np.sqrt(r * gRFT)      # km/s

    # Also plot pure Newtonian baryons for comparison
    v_Newton = np.sqrt(r * gN)

    ax.plot(r, v_RFT, label="RFT (baryons only)")
    ax.plot(r, v_Newton, linestyle="--", alpha=0.5, label="Newtonian baryons")

    # Highlight the 10–40 kpc band and ~220–240 km/s region
    ax.axvspan(10.0, 40.0, color="grey", alpha=0.1, label="10–40 kpc")
    ax.axhspan(220.0, 240.0, color="orange", alpha=0.1, label="220–240 km/s")

    ax.set_xlabel(r"Radius $r$ (kpc)")
    ax.set_ylabel(r"$v_{\rm rot}$ (km/s)")
    ax.set_title("RFT Rotation Curve (Γ_eff–linked a$_0$, baryons only)")
    ax.set_xlim(0, 60)
    ax.set_ylim(0, max(280, v_RFT.max() * 1.05))
    ax.grid(True)
    ax.legend()

    return ax


# ------------------------------------------------------------
# MACROSCOPIC COHERENCE / COLLAPSE (Ξ(t))
# ------------------------------------------------------------

def lambda_rft(Gamma_eff: float,
               m_g: float = 1.0,
               dx_um: float = 1.0) -> float:
    """
    RFT collapse rate:

      λ_RFT = λ0_ref * (Gamma_eff/GAMMA_REF) * (m/m_ref)^2 * (dx/dx_ref)^2

    with m_ref = 1 g, dx_ref = 1 μm.
    """
    factor_m = (m_g / M_REF_G) ** 2
    factor_x = (dx_um / DX_REF_UM) ** 2
    lam = LAMBDA0_REF * (Gamma_eff / GAMMA_REF) * factor_m * factor_x
    return lam


def coherence_plot(Gamma_eff: float, ax=None):
    """
    Plot Ξ(t) = 1 - exp(-λ t) for a 1 g, 1 mm, 1 μm silicon crystal.

    We show time from ~1 s up to 1e9 s on a log scale.
    """
    if ax is None:
        fig, ax = plt.subplots()

    # Times from 1 s to 1e9 s (log spaced)
    t = np.logspace(0, 9, 300)  # seconds

    lam = lambda_rft(Gamma_eff, m_g=1.0, dx_um=1.0)
    Xi = 1.0 - np.exp(-lam * t)

    ax.semilogx(t, Xi, label=fr"Ξ(t), λ≈{lam:.2e} s$^{{-1}}$")

    # Mark t = 1e7 s (experimental bound scale)
    t_bound = 1.0e7
    Xi_bound = 1.0 - math.exp(-lam * t_bound)
    ax.axvline(t_bound, linestyle="--", alpha=0.4,
               label=fr"$t=10^7$ s, Ξ≈{Xi_bound:.2e}")

    ax.set_xlabel(r"Time $t$ (s)")
    ax.set_ylabel(r"Ξ(t)")
    ax.set_title("RFT Macroscopic Coherence (1 g, 1 μm, Γ_eff–scaled λ)")
    ax.set_ylim(0, 1.0)
    ax.grid(True, which="both", ls=":")
    ax.legend()

    return ax


# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Visualize RFT outputs (CMB peaks, BAO, rotation curves, "
            "and coherence) using an evolved Γ_eff value."
        )
    )
    parser.add_argument(
        "--input",
        default="rft_evolution.npz",
        help="Path to evolution NPZ file containing Gamma_eff (default: rft_evolution.npz)",
    )
    parser.add_argument(
        "--save",
        metavar="PATH",
        help="Optional path to save the rendered figure instead of only showing it.",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Generate plots without opening an interactive window.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    Gamma_eff = load_gamma_eff(args.input)

    print(f"[plot_results] Γ_eff (from evolution) = {Gamma_eff:.6e}")

    r_bao = bao_scale_mpc(Gamma_eff)
    print(f"[plot_results] r_BAO(Γ_eff) ≈ {r_bao:.3f} Mpc")

    lam_lab = lambda_rft(Gamma_eff, m_g=1.0, dx_um=1.0)
    print(f"[plot_results] λ_RFT (1 g, 1 μm) ≈ {lam_lab:.3e} s^-1")

    # Set up 2x2 figure: CMB, P(k), rotation curve, coherence
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    # Top-left: CMB TT
    cmb_tt_template(Gamma_eff, ax=axes[0, 0])

    # Top-right: P(k) with BAO
    bao_plot(Gamma_eff, ax=axes[0, 1])

    # Bottom-left: rotation curve
    rotation_curve_rft(Gamma_eff, ax=axes[1, 0])

    # Bottom-right: coherence Ξ(t)
    coherence_plot(Gamma_eff, ax=axes[1, 1])

    fig.tight_layout()

    if args.save:
        out_path = pathlib.Path(args.save)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=300)
        print(f"[plot_results] Saved figure to {out_path.resolve()}")

    if not args.no_show:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
