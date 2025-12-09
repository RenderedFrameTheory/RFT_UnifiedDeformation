"""
RFT_UnifiedDeformation package.

This package implements the Rendered Frame Theory (RFT) unified deformation
cosmology model: background evolution of a five–field scalar system,
derivation of a single deformation scalar Gamma_eff, and its use in the
calculation of BAO, CMB acoustic templates, galaxy rotation curves and
macroscopic collapse rates.
"""

from .parameters import (
    RFTParameters,
    DEFAULT_PARAMS,
)

from .cosmology_model import (
    rft_background_rhs,
    initial_background_state,
    compute_energy_densities,
    compute_gamma_eff,
)

from .evolution import (
    run_background_evolution,
    save_background_solution,
    load_background_solution,
)

from .observables import (
    bao_scale_mpc,
    acoustic_peaks_ell,
    rotation_curve_rft,
    lambda_rft,
    Xi_of_t,
)

from .plots import (
    plot_cmb_tt_template,
    plot_pk_with_bao,
    plot_rotation_curve,
    plot_coherence_curve,
)

from .comparison_lcdm import (
    lcdm_rs_mpc,
    lcdm_ell_A,
    lcdm_rotation_curve_newtonian,
)
