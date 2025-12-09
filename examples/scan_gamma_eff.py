"""
Scan the RFT deformation sector and show how Gamma_eff, BAO scale,
galaxy acceleration scale a0 and collapse rate lambda_RFT co–move.

This script varies a single deformation–sensitive parameter and repeats
the background integration to illustrate the structural linkage in the
unified deformation model.
"""

import copy

import numpy as np

from rft import DEFAULT_PARAMS, run_background_evolution
from rft.observables import bao_scale_mpc, lambda_rft


def main():
    base_params = DEFAULT_PARAMS

    # Choose a parameter that affects the deformation field evolution.
    # Here we scale m_Gamma around its default value.
    mGamma_values = [0.7, 0.8, 0.9, 1.0, 1.1]

    print("m_Gamma  Gamma_eff        r_BAO[Mpc]     a0[(km/s)^2/kpc]   lambda_RFT(1g,1µm)[s^-1]")
    for mGamma in mGamma_values:
        params = copy.deepcopy(base_params)
        params.m_Gamma = mGamma

        t_arr, y_arr, Gamma_eff = run_background_evolution(params=params)
        r_bao = bao_scale_mpc(Gamma_eff, params)
        a0 = params.a0_ref * (Gamma_eff / params.Gamma_ref)
        lam = lambda_rft(Gamma_eff, m_g=params.m_ref_g, dx_um=params.dx_ref_um, params=params)

        print(f"{mGamma:7.3f}  {Gamma_eff:10.4e}  {r_bao:12.3f}  {a0:18.3f}  {lam:22.3e}")


if __name__ == "__main__":
    main()
