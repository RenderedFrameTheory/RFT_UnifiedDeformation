# RFT Cosmology v1 — Unified Deformation Model

This document describes the scalar–tensor cosmology implemented in this repository. The model is complete at the background level: field content, action ingredients, FRW equations and the definition and role of the deformation scalar \(\Gamma_{\text{eff}}\) are all specified.

## 1. Field content and metrics

Spacetime metric:

- Einstein–frame metric \(g_{\mu\nu}\) with signature \((-+++)\),
- spatially flat FRW background:
  \[
  ds^2 = -dt^2 + a^2(t)\,d\vec x^2.
  \]

Scalar fields:

- Five homogeneous scalar fields
  \[
  \varphi^I = \{\Phi,\Gamma,R,\Xi,\Psi\}, \quad I = 1,\dots,5.
  \]

Matter coupling:

- Matter and radiation couple to a conformally related metric
  \[
  \tilde g_{\mu\nu} = A^2(\Phi,R)\,g_{\mu\nu},
  \]
  with
  \[
  A(\Phi,R) = 1 + a_\Phi \Phi + a_R R.
  \]

Effective Planck mass:

- The gravitational coupling is modified by the deformation field \(\Gamma\) through
  \[
  M_*^2(\Gamma) = M_{\rm Pl}^2\,(1 + \xi_\Gamma \Gamma).
  \]

In the numerical implementation, overall powers of \(M_{\rm Pl}\) are absorbed into the choice of units, so only the dimensionless factor \(1 + \xi_\Gamma\Gamma\) appears.

## 2. Scalar potential

The scalar sector is governed by a quadratic potential with three cross–couplings:

\[
V(\Phi,\Gamma,R,\Xi,\Psi) =
\frac{1}{2} m_\Phi^2 \Phi^2
+ \frac{1}{2} m_\Gamma^2 \Gamma^2
+ \frac{1}{2} m_R^2 R^2
+ \frac{1}{2} m_\Xi^2 \Xi^2
+ \frac{1}{2} m_\Psi^2 \Psi^2
+ \lambda_{\Gamma R} \Gamma R
+ \lambda_{\Phi\Gamma} \Phi\Gamma
+ \lambda_{\Xi\Psi} \Xi\Psi
+ \Lambda_{\rm RFT}.
\]

All masses \(m_I\), couplings \(\lambda_{IJ}\) and the constant term \(\Lambda_{\rm RFT}\) are specified numerically in `rft/parameters.py`.

## 3. Action and matter sector

The full action can be written schematically as

\[
S = \int d^4x \sqrt{-g}\left[
\frac{1}{2} M_*^2(\Gamma) R(g)
- \frac{1}{2} \sum_{I=1}^5 g^{\mu\nu} \partial_\mu \varphi^I \partial_\nu \varphi^I
- V(\Phi,\Gamma,R,\Xi,\Psi)
\right]
+ S_m[\tilde g_{\mu\nu}, \psi_m],
\]

where:

- \(R(g)\) is the Ricci scalar of \(g_{\mu\nu}\),
- \(S_m\) is the matter and radiation action, minimally coupled to \(\tilde g_{\mu\nu} = A^2(\Phi,R) g_{\mu\nu}\),
- \(\psi_m\) denotes the matter fields.

The resulting background equations are implemented explicitly in `rft/cosmology_model.py`.

## 4. Background FRW equations

### 4.1 State vector

The background evolution is encoded in the state vector

\[
y = \big(a, H,
 \Phi, \dot\Phi,
 \Gamma, \dot\Gamma,
 R, \dot R,
 \Xi, \dot\Xi,
 \Psi, \dot\Psi\big).
\]

The ODE system evolves \(y(t)\) from an early initial time \(t=0\) to a later time \(t=t_{\rm end}\) chosen for numerical convenience.

### 4.2 Energy densities

With the conformal factor \(A(\Phi,R)\), the Einstein–frame energy densities of baryons and radiation are

\[
\rho_b(a,\Phi,R) = A^4(\Phi,R) \frac{\rho_{b0}}{a^3},
\quad
\rho_\gamma(a,\Phi,R) = A^4(\Phi,R) \frac{\rho_{\gamma 0}}{a^4}.
\]

The scalar–sector energy density is

\[
\rho_{\rm RFT} = K_\varphi + V,
\quad
K_\varphi = \frac{1}{2}\left(
\dot\Phi^2 + \dot\Gamma^2 + \dot R^2 + \dot\Xi^2 + \dot\Psi^2
\right).
\]

The total Einstein–frame energy density is

\[
\rho_{\rm tot} = \rho_b + \rho_\gamma + \rho_{\rm RFT}.
\]

These expressions are used in `energy_densities()` in `cosmology_model.py`.

### 4.3 Effective Planck mass

The effective Planck mass squared is

\[
M_*^2(\Gamma) = 1 + \xi_\Gamma \Gamma
\]

in the chosen units. Its time derivative is

\[
\dot M_*^2 = \xi_\Gamma \dot\Gamma.
\]

### 4.4 Friedmann and Raychaudhuri equations (approximate form)

The full covariant equations include terms from \(\dot M_*^2\) and \(\ddot M_*^2\). In this minimal implementation the background Hubble rate is evolved using an approximate Raychaudhuri equation that captures the main dependence on the matter and scalar sector:

\[
-2 M_*^2(\Gamma)\,\dot H \approx \rho_b + \frac{4}{3}\rho_\gamma + 2 K_\varphi.
\]

This is the expression coded in `rft_background_rhs()`. The combination of this equation with the Hubble definition \(\dot a = aH\) and the scalar equations gives a closed system for \(a(t)\) and the fields.

### 4.5 Scalar equations of motion

The scalar fields satisfy

\[
\ddot\Phi + 3H\dot\Phi + \frac{\partial V}{\partial\Phi} = S_\Phi,
\]
\[
\ddot\Gamma + 3H\dot\Gamma + \frac{\partial V}{\partial\Gamma}
- \frac{1}{2}\xi_\Gamma R(g) = S_\Gamma,
\]
\[
\ddot R + 3H\dot R + \frac{\partial V}{\partial R} = S_R,
\]
\[
\ddot\Xi + 3H\dot\Xi + \frac{\partial V}{\partial\Xi} = 0,
\]
\[
\ddot\Psi + 3H\dot\Psi + \frac{\partial V}{\partial\Psi} = 0.
\]

Here:

- \(R(g) = 6(2H^2 + \dot H)\) is the Ricci scalar of the FRW metric,
- the matter sources for \(\Phi\) and \(R\) are
  \[
  S_\Phi = -a_\Phi\,T_m,
  \quad
  S_R = -a_R\,T_m,
  \]
  with \(T_m \approx -\rho_b\) for pressureless baryons,
- in this implementation \(S_\Gamma = 0\) at the background level.

The derivatives of the potential are

\[
\frac{\partial V}{\partial\Phi} = m_\Phi^2 \Phi + \lambda_{\Phi\Gamma}\Gamma,
\]
\[
\frac{\partial V}{\partial\Gamma} = m_\Gamma^2 \Gamma + \lambda_{\Gamma R} R + \lambda_{\Phi\Gamma}\Phi,
\]
\[
\frac{\partial V}{\partial R} = m_R^2 R + \lambda_{\Gamma R}\Gamma,
\]
\[
\frac{\partial V}{\partial\Xi} = m_\Xi^2 \Xi + \lambda_{\Xi\Psi}\Psi,
\]
\[
\frac{\partial V}{\partial\Psi} = m_\Psi^2 \Psi + \lambda_{\Xi\Psi}\Xi.
\]

All of these expressions are implemented explicitly in `cosmology_model.py`.

## 5. Deformation scalar \(\Gamma_{\text{eff}}\)

The deformation field \(\Gamma(t)\) controls both the effective Planck mass \(M_*^2(\Gamma)\) and, through the background expansion, the sound horizon and distances that enter BAO and CMB observables.

To link this to a single scalar parameter that can be compared across models and used in simplified templates, the effective deformation scalar is defined as

\[
\Gamma_{\text{eff}} =
\frac{
\displaystyle \int_{\eta_i}^{\eta_f}
\Gamma(\eta)\,W(\eta)\,d\eta
}{
\displaystyle \int_{\eta_i}^{\eta_f}
W(\eta)\,d\eta
},
\quad
W(\eta) \propto a^2(\eta)\rho_b(\eta).
\]

In terms of cosmic time \(t\), using \(d\eta = dt/a\), this becomes

\[
\Gamma_{\text{eff}} =
\frac{
\displaystyle \int_{t_i}^{t_f}
\Gamma(t)\,a(t)\,\rho_b(t)\,dt
}{
\displaystyle \int_{t_i}^{t_f}
a(t)\,\rho_b(t)\,dt
}.
\]

The weight \(a\rho_b\) gives more importance to epochs when the baryon fluid dominates the acoustic and structure signals. This definition is implemented in `compute_gamma_eff()`.

## 6. Observables derived from \(\Gamma_{\text{eff}}\)

### 6.1 BAO scale

The BAO scale in real space is modelled as

\[
r_{\rm BAO}(\Gamma_{\text{eff}})
= C_{\rm BAO}\,\frac{\pi}{\sqrt{\Gamma_{\text{eff}}}},
\]

with the constant \(C_{\rm BAO}\) fixed by requiring

\[
r_{\rm BAO}(\Gamma_{\rm ref}) = r_{\rm BAO,ref},
\]

where \(\Gamma_{\rm ref}\) and \(r_{\rm BAO,ref}\) are given in the parameter set. This choice captures the inverse square–root dependence of the effective sound horizon on the deformation strength.

### 6.2 CMB acoustic peaks

The CMB acoustic peak positions are encoded in a template form:

- Reference peak multipoles at \(\Gamma_{\rm ref}\) are
  \[
  \ell_{\rm ref} = [220, 540, 800, 1100].
  \]
- For general \(\Gamma_{\text{eff}}\),
  \[
  \ell_n(\Gamma_{\text{eff}}) =
  \ell_{n,{\rm ref}} \sqrt{\frac{\Gamma_{\text{eff}}}{\Gamma_{\rm ref}}}.
  \]

The full CMB spectrum is not computed here. Instead, these peaks are used to construct a TT template made from Gaussians centred at \(\ell_n(\Gamma_{\text{eff}})\). This keeps the focus on how the acoustic scale shifts when \(\Gamma_{\text{eff}}\) changes.

### 6.3 Galaxy rotation curve

In the weak–field, quasi–static limit around a galaxy, the RFT model produces an effective acceleration law that can be approximated by a MOND–like simple–\(\mu\) form:

\[
g_{\rm RFT}(r) \approx
\frac{1}{2} \left[
g_N(r) + \sqrt{g_N^2(r) + 4 g_N(r) a_0(\Gamma_{\text{eff}})}
\right],
\]

where:

- \(g_N(r) = G M(<r)/r^2\) is the Newtonian acceleration from the baryons alone,
- \(a_0(\Gamma_{\text{eff}})\) is the RFT acceleration scale,
  \[
  a_0(\Gamma_{\text{eff}}) = a_{0,{\rm ref}} \frac{\Gamma_{\text{eff}}}{\Gamma_{\rm ref}}.
  \]

The circular velocity is

\[
v_{\rm RFT}^2(r) = r\,g_{\rm RFT}(r).
\]

With an exponential disk
\[
M(<r) = M_{\rm disk}\left[1 - \left(1 + \frac{r}{R_d}\right)e^{-r/R_d}\right],
\]
and an appropriate choice of \(a_{0,{\rm ref}}\), the RFT rotation curve flattens near \(220\)–\(240\) km/s from about 10 to 40 kpc, using baryons only and a single \(\Gamma_{\text{eff}}\).

### 6.4 Macroscopic collapse rate

The RFT macroscopic collapse rate is modelled as

\[
\lambda_{\rm RFT}(m,\Delta x;\Gamma_{\text{eff}}) =
\lambda_0 \left(\frac{\Gamma_{\text{eff}}}{\Gamma_{\rm ref}}\right)
\left(\frac{m}{m_0}\right)^2
\left(\frac{\Delta x}{\Delta x_0}\right)^2,
\]

with reference scales \(m_0, \Delta x_0, \lambda_0\) fixed in `parameters.py`.

The coherence–loss probability for a superposition lasting time \(t\) is

\[
\Xi(t) = 1 - e^{-\lambda_{\rm RFT} t}.
\]

For \(m = 1~{\rm g}\), \(\Delta x = 1~\mu{\rm m}\) and a reference value \(\lambda_0 \approx 10^{-10}~{\rm s}^{-1}\), \(\Xi(t)\) stays well below 1 for \(t \lesssim 10^7~{\rm s}\), in line with current macroscopic coherence bounds, while still allowing the collapse rate to grow if \(\Gamma_{\text{eff}}\) differs from \(\Gamma_{\rm ref}\).

## 7. Unified role of \(\Gamma_{\text{eff}}\)

The key structural feature of this model is that \(\Gamma_{\text{eff}}\), derived once from the background evolution, enters all sectors in a fixed way:

- BAO:
  \[
  r_{\rm BAO}(\Gamma_{\text{eff}}) \propto \frac{1}{\sqrt{\Gamma_{\text{eff}}}}.
  \]
- CMB acoustic scale:
  \[
  \ell_n(\Gamma_{\text{eff}}) \propto \sqrt{\Gamma_{\text{eff}}}.
  \]
- Galaxy rotation:
  \[
  a_0(\Gamma_{\text{eff}}) \propto \Gamma_{\text{eff}}.
  \]
- Collapse:
  \[
  \lambda_{\rm RFT} \propto \Gamma_{\text{eff}} m^2 \Delta x^2.
  \]

Adjusting the parameters that control the evolution of \(\Gamma(t)\) to benefit cosmology automatically changes the acceleration scale in rotation curves and the strength of macroscopic collapse. There is no independent tuning of BAO, galaxy dynamics and collapse; they share one deformation scalar produced by one cosmological core.
