import argparse
import pathlib

import matplotlib.pyplot as plt
import numpy as np

from rft import (
    DEFAULT_PARAMS,
    load_background_solution,
)
from rft.observables import (
    bao_scale_mpc,
    acoustic_peaks_ell,
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
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Visualise RFT unified deformation outputs using a saved "
            "background NPZ file containing Gamma_eff."
        )
    )
    parser.add_argument(
        "--input",
        default="rft_background.npz",
        help="Path to NPZ file with background evolution and Gamma_eff "
             "(default: rft_background.npz)",
    )
    parser.add_argument(
        "--save",
        metavar="PATH",
        help="Optional path to save the 2x2 diagnostics figure.",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not open an interactive window; useful for batch runs.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    npz_path = pathlib.Path(args.input)
    if not npz_path.exists():
        raise FileNotFoundError(
            f"Input NPZ file not found at '{npz_path}'. "
            "Run examples/run_full_pipeline.py first or specify --input."
        )

    # Load background and Gamma_eff
    t_arr, y_arr, Gamma_eff, extras = load_background_solution(str(npz_path))

    params = DEFAULT_PARAMS

    print(f"[RFT] Loaded background from: {npz_path.resolve()}")
    print(f"[RFT] Gamma_eff = {Gamma_eff:.6e}")

    # Core RFT observables from Gamma_eff
    r_bao = bao_scale_mpc(Gamma_eff, params)
    ell_peaks = acoustic_peaks_ell(Gamma_eff, params)
    a0 = params.a0_ref * (Gamma_eff / params.Gamma_ref)
    lam_ref = lambda_rft(
        Gamma_eff,
        m_g=params.m_ref_g,
        dx_um=params.dx_ref_um,
        params=params,
    )

    print(f"[RFT] r_BAO(Gamma_eff) ≈ {r_bao:.3f} Mpc")
    print(f"[RFT] Acoustic peak template ell_n(Gamma_eff) = {ell_peaks}")
    print(f"[RFT] a0(Gamma_eff) ≈ {a0:.3f} (km/s)^2/kpc")
    print(f"[RFT] lambda_RFT(1 g, 1 µm) ≈ {lam_ref:.3e} s^-1")

    # Simple ΛCDM reference
    rs_lcdm = lcdm_rs_mpc()
    D_A_ref = 14000.0  # rough reference angular diameter distance in Mpc
    ellA_lcdm = lcdm_ell_A(D_A_ref, rs_lcdm)
    print(f"[ΛCDM] Approximate r_s ≈ {rs_lcdm:.3f} Mpc")
    print(f"[ΛCDM] Approximate ell_A ≈ {ellA_lcdm:.1f} (with D_A ≈ {D_A_ref:.1f} Mpc)")

    # Build 2x2 diagnostics figure
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    # Top-left: CMB TT template
    plot_cmb_tt_template(Gamma_eff, params=params, ax=axes[0, 0])

    # Top-right: P(k) with BAO feature
    plot_pk_with_bao(Gamma_eff, params=params, ax=axes[0, 1])

    # Bottom-left: rotation curve (RFT vs pure baryons)
    plot_rotation_curve(Gamma_eff, params=params, ax=axes[1, 0])

    # Bottom-right: macroscopic coherence Xi(t)
    plot_coherence_curve(
        Gamma_eff,
        m_g=params.m_ref_g,
        dx_um=params.dx_ref_um,
        params=params,
        ax=axes[1, 1],
    )

    fig.tight_layout()

    if args.save:
        out_path = pathlib.Path(args.save)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=300)
        print(f"[RFT] Saved diagnostics figure to {out_path.resolve()}")

    if not args.no_show:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
