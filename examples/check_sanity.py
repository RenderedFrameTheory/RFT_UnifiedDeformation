from rft import DEFAULT_PARAMS, load_background_solution
from rft.observables import bao_scale_mpc, acoustic_peaks_ell, lambda_rft

def main():
    t, y, Gamma_eff, extras = load_background_solution("rft_background.npz")
    p = DEFAULT_PARAMS

    r_bao = bao_scale_mpc(Gamma_eff, p)
    ells = acoustic_peaks_ell(Gamma_eff, p)
    lam  = lambda_rft(Gamma_eff, m_g=p.m_ref_g, dx_um=p.dx_ref_um, params=p)

    print(f"Gamma_eff = {Gamma_eff:.6e}")
    print(f"BAO scale ≈ {r_bao:.2f} Mpc")
    print(f"Peak ells ≈ {', '.join(f'{e:.0f}' for e in ells)}")
    print(f"lambda_RFT(1 g, 1 μm) ≈ {lam:.3e} s^-1")

if __name__ == '__main__':
    main()
