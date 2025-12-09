# RFT Unified Deformation — Predictions and Tests

This document summarises the main qualitative predictions of the RFT unified deformation model and the types of tests that can distinguish it from the standard ΛCDM+MOND+CSL stack.

## 1. Co–movement of BAO, galaxy acceleration and collapse

### 1.1 Linked shifts

The model predicts that:

- If the BAO scale \(r_{\rm BAO}\) is inferred to differ from the ΛCDM expectation by a factor \(f_{\rm BAO}\),
  \[
  r_{\rm BAO} \propto \frac{1}{\sqrt{\Gamma_{\text{eff}}}},
  \]
  then
  \[
  \Gamma_{\text{eff}} \propto \frac{1}{r_{\rm BAO}^2}.
  \]

- The galaxy acceleration scale behaves as
  \[
  a_0(\Gamma_{\text{eff}}) \propto \Gamma_{\text{eff}},
  \]
  so a relative change in \(r_{\rm BAO}\) should correspond to a linked change in \(a_0\).

- The macroscopic collapse rate scales as
  \[
  \lambda_{\rm RFT} \propto \Gamma_{\text{eff}} m^2 \Delta x^2.
  \]

Taken together, this means any deviation in the inferred BAO scale that is accommodated by adjusting \(\Gamma_{\text{eff}}\) must be accompanied by predictable shifts in galaxy rotation and collapse constraints.

### 1.2 Testable correlation

This implies a test:

- Measure BAO scale and CMB acoustic scale at the best available precision.
- Determine \(\Gamma_{\text{eff}}\) from the cosmological data within the RFT framework.
- Use this \(\Gamma_{\text{eff}}\) to compute a predicted \(a_0\).
- Compare the predicted acceleration scale with the value inferred from galaxy rotation curves.
- Check whether macroscopic coherence experiments are compatible with the corresponding \(\lambda_{\rm RFT}\).

A consistent RFT fit must pass all three checks with a single \(\Gamma_{\text{eff}}\).

## 2. Baryon–only rotation curves and acceleration scale

RFT predicts that Milky Way–like galaxies can be fitted with baryons only and a deformation–linked acceleration scale \(a_0(\Gamma_{\text{eff}})\), using the same \(\Gamma_{\text{eff}}\) that controls BAO and CMB.

A falsifiable test is:

- Fix RFT parameters using cosmology alone to determine \(\Gamma_{\text{eff}}\).
- Compute \(a_0(\Gamma_{\text{eff}})\) with no further tuning.
- Apply the RFT rotation law to a suite of galaxies with known baryon distributions.
- Check whether the rotation curves are reproduced as well as, or better than, standard MOND with a fitted a₀ or NFW halos.

If cosmology–fixed RFT cannot reproduce the observed rotation curves within reasonable errors, the model fails as a unified description.

## 3. Collapse bounds and future tests

RFT predicts a specific dependence of the macroscopic collapse rate on \(\Gamma_{\text{eff}}\), mass and separation:

\[
\lambda_{\rm RFT}(m,\Delta x;\Gamma_{\text{eff}}) =
\lambda_0 \left(\frac{\Gamma_{\text{eff}}}{\Gamma_{\rm ref}}\right)
\left(\frac{m}{m_0}\right)^2
\left(\frac{\Delta x}{\Delta x_0}\right)^2.
\]

Given \(\Gamma_{\text{eff}}\) from cosmology and \(\lambda_0\) fixed by one reference experiment, RFT makes definite predictions for:

- Larger masses,
- Larger separations,
- Longer coherence times.

A decisive test is to design experiments in a regime where:

- Standard CSL predictions differ from decoherence–only predictions,
- RFT's scaling in \(\Gamma_{\text{eff}}\) leads to a different dependence on mass and separation.

If experiments choose parameters where CSL and RFT diverge, a consistent set of results across several setups can support or exclude the RFT scaling.

## 4. Relation to ΛCDM tensions

The model is particularly relevant if current cosmological tensions persist or grow:

- If BAO or CMB–inferred distances deviate from simple ΛCDM fits in a way that suggests a modified sound horizon, a non–standard \(\Gamma_{\text{eff}}\) may be preferred.
- RFT then predicts a linked adjustment in galaxy acceleration and collapse scale, not just a local cosmology fix.

The code in this repository does not claim to resolve specific tensions numerically. It provides the structure required to make such links explicit once a full numerical fit is attempted.

## 5. Summary

The core prediction of this implementation is structural:

- There is one deformation scalar \(\Gamma_{\text{eff}}\),
- The BAO scale, CMB acoustic scale, galaxy acceleration and macroscopic collapse rate are all functions of this scalar,
- Adjusting the scalar to match one class of observations imposes specific, testable consequences in the others.

The numerical values in this repository are not tuned to data. They are chosen to give realistic orders of magnitude and stable evolution. The predicted correlations are, however, fixed by the formulas in the code and can be confronted with data once a full numerical analysis is implemented.
