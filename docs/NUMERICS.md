# Numerical Methods

This document explains the numerical scheme used in the RFT_UnifiedDeformation code.

## 1. Background evolution

The background system is a set of 12 first–order ODEs:

- Variables:
  \[
  y = (a, H,
       \Phi, \dot\Phi,
       \Gamma, \dot\Gamma,
       R, \dot R,
       \Xi, \dot\Xi,
       \Psi, \dot\Psi).
  \]

- Right–hand side: `rft_background_rhs(t, y, params)` in `cosmology_model.py`.

The equations are integrated with SciPy's `solve_ivp`:

- Method: `RK45` (explicit Dormand–Prince 4(5) method).
- Relative tolerance: `rtol = 1e-7`.
- Absolute tolerance: `atol = 1e-9`.
- Time interval: `t_start = 0`, `t_end = 1e4` in core units.
- Evaluation grid: `n_steps` equally spaced points in `[t_start, t_end]`.

The initial state is constructed from the parameters in `initial_background_state()`.

## 2. Calculation of \(\Gamma_{\text{eff}}\)

After the background integration:

- The code extracts arrays:
  - `t_arr` (time),
  - `Gamma_arr` (deformation field),
  - `a_arr` (scale factor),
  - `rho_b_arr` (baryon density from `compute_energy_densities`).

- \(\Gamma_{\text{eff}}\) is computed via `compute_gamma_eff()` using the trapezoidal rule:

  - Weight \(W(t) = a(t)\rho_b(t)\),
  - Numerator: \(\sum \Gamma_{\rm mid} W_{\rm mid} \Delta t\),
  - Denominator: \(\sum W_{\rm mid} \Delta t\),
  - Midpoint values are simple averages between adjacent grid points.

The result is returned as a scalar `Gamma_eff`.

## 3. Saving and loading background data

The background solution and derived quantities are saved to a `.npz` file using `save_background_solution()`:

- Keys:
  - `t`, `y`, `Gamma_eff`,
  - `a`, `Phi`, `Gamma`, `R`, `Xi`, `Psi`,
  - `Phi_dot`, `Gamma_dot`, `R_dot`, `Xi_dot`, `Psi_dot`,
  - `rho_b`, `rho_gamma`, `rho_rft`, `rho_tot`.

These arrays can be reloaded with `load_background_solution(path)`.

## 4. Observable calculations

Observable–level calculations use simple analytic or template forms:

- `bao_scale_mpc(Gamma_eff)`:
  - Computes the BAO scale from \(\Gamma_{\text{eff}}\) using a fixed inverse square–root law with one calibration constant.

- `acoustic_peaks_ell(Gamma_eff)`:
  - Returns four peak multipoles from a square–root scaling with \(\Gamma_{\text{eff}}\).

- `rotation_curve_rft(r_kpc, Gamma_eff)`:
  - Computes the Newtonian acceleration from an exponential baryon disk,
  - Computes the RFT acceleration using a MOND–like formula with \(a_0(\Gamma_{\text{eff}})\),
  - Returns RFT and Newtonian rotation velocities.

- `lambda_rft(Gamma_eff, m_g, dx_um)` and `Xi_of_t(t, Gamma_eff, m_g, dx_um)`:
  - Provide the macroscopic collapse rate and coherence–loss curve.

These functions use only algebraic expressions and do not involve additional numerical integration.

## 5. Syntax for running examples

The `examples` directory contains two scripts:

- `run_full_pipeline.py`:
  - Runs `run_background_evolution()`,
  - Saves the background solution,
  - Prints \(\Gamma_{\text{eff}}\) and the derived BAO scale and acceleration scale,
  - Generates four plots using functions in `rft/plots.py`.

- `scan_gamma_eff.py`:
  - Perturbs one of the parameters affecting the deformation field,
  - Repeats the background evolution for a small set of parameter values,
  - Shows how \(\Gamma_{\text{eff}}\), \(r_{\rm BAO}\), \(a_0\) and \(\lambda_{\rm RFT}\) move together.

Both scripts can be run directly with Python once the requirements are installed.

The numerical scheme is deliberately simple and explicit. It is sufficient to explore the structural link between the deformation scalar and the four classes of observables. It is not a substitute for a full precision Boltzmann solver or galaxy rotation curve catalog fit.
