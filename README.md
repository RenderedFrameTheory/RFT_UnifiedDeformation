# RFT_UnifiedDeformation

Rendered Frame Theory (RFT) Cosmology — Unified Deformation Model

This repository implements a single, explicit scalar–tensor cosmology where one deformation scalar \(\Gamma_{\text{eff}}\) is derived from a five–field core and then used to tie together four sectors:

- BAO sound horizon / BAO scale
- CMB acoustic scale and peak positions (template level)
- Baryon–only galaxy rotation curves (no dark matter halo)
- Macroscopic collapse rate for laboratory superpositions

There is one RFT deformation engine, one derived \(\Gamma_{\text{eff}}\), and four classes of observables forced to move together. There is no separate tuning of cosmology, rotation curves, and collapse.

The model is written in standard scalar–tensor language:

- A FRW background with variable effective Planck mass \(M_*^2(\Gamma)\).
- Five homogeneous scalar fields \(\Phi,\Gamma,R,\Xi,\Psi\) with a fixed quadratic–coupling potential.
- Matter coupled conformally to a physical metric \(\tilde g_{\mu\nu} = A^2(\Phi,R) g_{\mu\nu}\).
- A deformation field \(\Gamma\) that controls both gravity strength and the cosmological sound horizon.

The code provides:

- A background evolution solver for the five–field RFT system.
- A definition and numerical evaluation of \(\Gamma_{\text{eff}}\).
- Functions to compute BAO scale, CMB acoustic template peaks, a MOND–like RFT galaxy rotation curve, and a macroscopic collapse curve from \(\Gamma_{\text{eff}}\).
- A minimal ΛCDM comparison module for honest side–by–side checks.

This is a self–contained, inspectable testbed. It is not a full Boltzmann or N–body pipeline. The approximations used are explicit and documented. The aim is structural unification: to show that one deformation scalar, produced by one core, can consistently feed four observational sectors without per–sector retuning.

## Repository layout

```text
RFT_UnifiedDeformation/
├── README.md
├── requirements.txt
├── .gitignore
├── rft/
│   ├── __init__.py
│   ├── parameters.py
│   ├── cosmology_model.py
│   ├── evolution.py
│   ├── observables.py
│   ├── plots.py
│   └── comparison_lcdm.py
├── docs/
│   ├── THEORY.md
│   ├── PARAMETERS.md
│   ├── MODEL_COMPARISON.md
│   ├── PREDICTIONS.md
│   └── NUMERICS.md
└── examples/
    ├── run_full_pipeline.py
    └── scan_gamma_eff.py
Installation

Create and activate a virtual environment, then install the dependencies:

pip install -r requirements.txt

Quick start

Run a full background evolution, compute 
Γ
eff
Γ
eff
	​

, and generate all four plots plus a ΛCDM comparison:

python examples/run_full_pipeline.py


Scan over small deformations of the parameters to see how 
Γ
eff
Γ
eff
	​

, the BAO scale, the galaxy acceleration scale 
𝑎
0
a
0
	​

, and the collapse rate 
𝜆
RFT
λ
RFT
	​

 co–move:

python examples/scan_gamma_eff.py


All parameters used, and their roles, are documented in docs/PARAMETERS.md. The theoretical structure is detailed in docs/THEORY.md. A direct comparison with ΛCDM, MOND–like gravity and CSL–style collapse models is given in docs/MODEL_COMPARISON.md.
## What this *is* and what this *is not*

- This *is* a minimal, inspectable RFT cosmology core where a single Γ_eff
  drives BAO scale, CMB acoustic peaks, baryon-only rotation curves and
  macroscopic collapse in one shot.

- This is *not*:
  - a full Boltzmann solver,
  - a parameter-inference engine,
  - a claim of beating ΛCDM on Planck χ².

If you want to attack it, attack the actual equations in `rft/cosmology_model.py`
and the Γ_eff → {BAO, ℓ, a₀, λ} mappings in `rft/observables.py`.
