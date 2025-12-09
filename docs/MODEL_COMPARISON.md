# RFT Unified Deformation vs ΛCDM, MOND and CSL

This document compares the structure of the RFT unified deformation model with the usual combination of ΛCDM, MOND–like gravity and CSL–style collapse models.

## 1. Components of the standard stack

The standard picture usually combines:

1. **ΛCDM cosmology**

   - Parameters: \(\Omega_b, \Omega_c, \Omega_\Lambda, H_0, n_s, A_s\).
   - Governs background expansion, CMB, BAO, matter power spectrum.

2. **Galaxy–scale gravity**

   - Either dark matter halos (NFW profiles) with masses and concentrations fitted per galaxy,
   - or MOND–like interpolation with a universal acceleration \(a_0\).

3. **Collapse and decoherence**

   - Often pure environmental decoherence.
   - In objective collapse proposals, additional parameters (\(\lambda, r_c\) for CSL).

These components are logically and parametrically separate. Changing, for example, the MOND acceleration \(a_0\) does not have a prescribed effect on the BAO scale or CMB peaks. Adjusting CSL parameters does not feed back into ΛCDM expansion.

## 2. Components of the RFT unified deformation model

The RFT model implemented here consists of:

1. **Single scalar–tensor cosmology**

   - Five scalar fields \(\Phi,\Gamma,R,\Xi,\Psi\),
   - Conformal matter coupling via \(A(\Phi,R)\),
   - Effective Planck mass via \(M_*^2(\Gamma)\),
   - Quadratic potential with fixed couplings,
   - Background FRW evolution with all fields.

2. **Single deformation scalar**

   - The deformation field \(\Gamma(t)\) is evolved along with the other fields.
   - A single scalar \(\Gamma_{\text{eff}}\) is defined as a weighted average over the background:
     \[
     \Gamma_{\text{eff}} =
     \frac{\int \Gamma(t) a(t) \rho_b(t) dt}{\int a(t) \rho_b(t) dt}.
     \]

3. **Four sectors tied to \(\Gamma_{\text{eff}}\)**

   - BAO scale \(r_{\rm BAO}(\Gamma_{\text{eff}})\),
   - CMB acoustic template peaks \(\ell_n(\Gamma_{\text{eff}})\),
   - Galaxy acceleration scale \(a_0(\Gamma_{\text{eff}})\),
   - Macroscopic collapse rate \(\lambda_{\rm RFT}(\Gamma_{\text{eff}},m,\Delta x)\).

There are no separate deformation parameters for each sector.

## 3. Parameter structure

### Standard stack

- ΛCDM: ~6 global parameters.
- MOND or halo profiles:
  - MOND: one acceleration \(a_0\) plus a function choice \(\mu(x)\),
  - Halos: several parameters per galaxy (mass, concentration, etc.).
- CSL collapse:
  - At least two new parameters \(\lambda, r_c\).

The stack uses a comparable total number of parameters as RFT but they are split across independent frameworks.

### RFT unified deformation

RFT uses:

- Scalar sector:
  - 5 masses \(m_I\),
  - 3 cross–couplings \(\lambda_{IJ}\),
  - 1 constant term \(\Lambda_{\rm RFT}\),
  - 3 metric couplings \(a_\Phi, a_R, \xi_\Gamma\).
- Background densities and initial conditions:
  - 2 density scales \(\rho_{b0}, \rho_{\gamma0}\),
  - initial field values and velocities,
  - initial \(a_0, H_0\).
- Calibration constants:
  - \(r_{\rm BAO,ref}\), \(a_{0,{\rm ref}}\), \(\lambda_0\).

The calibration constants are fixed at single reference points and then held fixed. They set units and absolute scales, not new functional degrees of freedom.

The main difference is that once the scalar sector parameters and initial conditions are chosen to produce a particular \(\Gamma_{\text{eff}}\), the same \(\Gamma_{\text{eff}}\) feeds all observables.

## 4. BAO and CMB

In ΛCDM:

- The BAO sound horizon and CMB peaks are determined by the matter–radiation content and expansion history encoded in \(\Omega_b, \Omega_c, \Omega_\Lambda, H_0\).
- Galaxy–scale parameters and collapse parameters do not enter those calculations.

In RFT:

- The sound horizon and acoustic scale depend on the same background expansion, but the deformation field \(\Gamma\) modifies the effective Planck mass and the detailed Hubble history.
- The effect of \(\Gamma\) over time is summarised by \(\Gamma_{\text{eff}}\), which then sets
  \[
  r_{\rm BAO}(\Gamma_{\text{eff}}) \propto \frac{1}{\sqrt{\Gamma_{\text{eff}}}},
  \quad
  \ell_n(\Gamma_{\text{eff}}) \propto \sqrt{\Gamma_{\text{eff}}}.
  \]

The code implements this relationship explicitly. It does not yet replace a full Boltzmann solver; instead, it provides a transparent mapping between \(\Gamma_{\text{eff}}\) and the main observational scales.

## 5. Galaxy rotation curves

In ΛCDM:

- Flat rotation curves are explained by dark matter halos,
- or, in MOND, by modifying the acceleration law with a universal \(a_0\).

The chosen value of \(a_0\) or halo parameters does not have a defined impact on BAO or CMB.

In RFT:

- The baryon–only rotation curve is computed from an exponential disk mass profile and an effective acceleration law
  \[
  g_{\rm RFT}(r) =
  \frac{1}{2}\left[g_N(r) + \sqrt{g_N^2(r) + 4 g_N(r) a_0(\Gamma_{\text{eff}})}\right],
  \]
  with
  \[
  a_0(\Gamma_{\text{eff}}) = a_{0,{\rm ref}} \frac{\Gamma_{\text{eff}}}{\Gamma_{\rm ref}}.
  \]
- Once \(\Gamma_{\text{eff}}\) and the baryon distribution are fixed, the flat part of the rotation curve is determined.

The same \(\Gamma_{\text{eff}}\) that sets \(r_{\rm BAO}\) sets \(a_0\). There is no independently tunable a₀.

## 6. Collapse and decoherence

In the standard stack:

- Decoherence is often handled by environmental models,
- Objective collapse, if used, introduces new parameters \(\lambda, r_c\) that are not linked to ΛCDM parameters.

In RFT:

- The collapse rate depends on \(\Gamma_{\text{eff}}\) and the mass and separation of the superposition:
  \[
  \lambda_{\rm RFT}(m,\Delta x;\Gamma_{\text{eff}}) \propto
  \Gamma_{\text{eff}} m^2 \Delta x^2.
  \]
- Changing \(\Gamma_{\text{eff}}\) to adjust BAO or rotation curves automatically changes \(\lambda_{\rm RFT}\).

This ties macroscopic coherence tests directly to the same deformation scalar that controls cosmology and galaxy dynamics.

## 7. Summary

RFT does not reduce the parameter count below the combined ΛCDM+MOND+CSL stack; that is not the claim. The difference is structural:

- ΛCDM+MOND+CSL uses several frameworks with separate parameter spaces.
- RFT uses one scalar–tensor core and one derived deformation scalar \(\Gamma_{\text{eff}}\) to control BAO, CMB acoustic scales, galaxy acceleration and collapse strength.

The code in this repository implements that structure explicitly. Any change in the scalar sector that affects \(\Gamma_{\text{eff}}\) is reflected, by design, in all four observable sectors. There is no hidden per–sector retuning.
