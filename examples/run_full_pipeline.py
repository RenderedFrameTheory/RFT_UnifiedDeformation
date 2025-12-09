"""
Run the full RFT unified deformation pipeline:

  - integrate the background evolution,
  - compute Gamma_eff,
  - compute BAO scale, acoustic peaks, a0 and collapse rate,
  - generate four diagnostic plots,
  - compare with simple ΛCDM reference scales.

This script is intended as the main entry point for inspecting the model.
"""

import pathlib

import numpy as np
import matplotlib.pyplot as plt

from rft import (
    DEFAULT_PARAMS,
    run_background_evolution,
    save_background_solution,
)
from rft.observables import (
    bao_scale_mpc,
    acoustic_peaks_ell,
    rotation_curve_rft,
    lambda_rft,
)
from rft.plots import (
    plot_cmb_tt_template,
    plot_pk_with_bao,
    plot_rotation_curve,
    plot_coherence_curve,
)
from rft.comparison_lcdm import (
    lcdm_rs_mpc,
    lcdm_ell_A,
    lcdm_rotation_curve_newtonian,
)


def main():
    params = DEFAULT_PARAMS

    # Integrate the background and compute Gamma_eff
    t_arr, y_arr, Gamma_eff = run_background_evolution(params=params)

    print(f"[RFT] Gamma_eff = {Gamma_eff:.6e}")

    # Save background solution
    out_path = pathlib.Path("rft_background.npz")
    save_background_solution(str(out_path), t_arr, y_arr, Gamma_eff, params=params)
    print(f"[RFT] Background solution saved to {out_path.resolve()}")

    # Derived observables
    r_bao = bao_scale_mpc(Gamma_eff, params)
    ell_peaks = acoustic_peaks_ell(Gamma_eff, params)
    a0 = params.a0_ref * (Gamma_eff / params.Gamma_ref)
    lam_ref = lambda_rft(Gamma_eff, m_g=params.m_ref_g, dx_um=params.dx_ref_um, params=params)

    print(f"[RFT] r_BAO(Gamma_eff) ≈ {r_bao:.3f} Mpc")
    print(f"[RFT] Acoustic peak template ell_n(Gamma_eff) = {ell_peaks}")
    print(f"[RFT] a0(Gamma_eff) ≈ {a0:.3f} (km/s)^2/kpc")
    print(f"[RFT] lambda_RFT(1 g, 1 µm) ≈ {lam_ref:.3e} s^-1")

    # Simple ΛCDM reference scales
    rs_lcdm = lcdm_rs_mpc()
    # Use a rough reference D_A for comparison (not integrated here)
    D_A_ref = 14000.0
    ellA_lcdm = lcdm_ell_A(D_A_ref, rs_lcdm)
    print(f"[ΛCDM] Approximate r_s ≈ {rs_lcdm:.3f} Mpc")
    print(f"[ΛCDM] Approximate ell_A ≈ {ellA_lcdm:.1f} (using D_A ≈ {D_A_ref:.1f} Mpc)")

    # Build plots
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    # CMB TT template
    plot_cmb_tt_template(Gamma_eff, params=params, ax=axes[0, 0])

    # P(k) with BAO
    plot_pk_with_bao(Gamma_eff, params=params, ax=axes[0, 1])

    # Rotation curve
    plot_rotation_curve(Gamma_eff, params=params, ax=axes[1, 0])

    # Coherence curve
    plot_coherence_curve(Gamma_eff, m_g=params.m_ref_g, dx_um=params.dx_ref_um,
                         params=params, ax=axes[1, 1])

    fig.tight_layout()
    fig.savefig("rft_unified_deformation_diagnostics.png", dpi=300)
    print("[RFT] Diagnostic plots saved to rft_unified_deformation_diagnostics.png")

    plt.show()


if __name__ == "__main__":
    main()
