**Author:** Dan Alec Yamaguchi (GitHub: danalec)
**Affiliation:** Independent Researcher
**Email:** danalec@gmail.com
**ORCID:** [0009-0002-9725-7779](https://orcid.org/0009-0002-9725-7779)
**DOI:** [10.5281/zenodo.20357668](https://doi.org/10.5281/zenodo.20357668)
**Date:** 22 May 2026

---

## Downloads

| File    | Link                                                         | SHA3-512                                                     |
| ------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| PDF     | [yamaguchi-rh-2026.pdf](./yamaguchi-rh-2026.pdf)             | `24f418634ba66e817c0e5c6e4a26d2e6cca0a450272b2ef644f79b82098c0ef9c6a9b3e5516f359c212f3381ca401352047ebe3643acce75a4fafbbd16896483` |
| TEX     | [yamaguchi-rh-2026.tex](./yamaguchi-rh-2026.tex)             | `c904eff4d06cd0b9503ccd176187dc01d2931e13fdf278487a6eac5e0b725f55371c39de80b80b75cddff30e118daf7e0e9eb7f52249683d593aa849c93dd4ec` |
| PDF.sig | [yamaguchi-rh-2026.pdf.sig](./yamaguchi-rh-2026.pdf.sig)     | Bitcoin signature                                            |
| TEX.sig | [yamaguchi-rh-2026.tex.sig](./yamaguchi-rh-2026.tex.sig)     | Bitcoin signature                                            |
| PDF.ots | [yamaguchi-rh-2026.pdf.ots](./yamaguchi-rh-2026.pdf.ots)     | Bitcoin timestamp                                            |
| TEX.ots | [yamaguchi-rh-2026.tex.ots](./yamaguchi-rh-2026.tex.ots)     | Bitcoin timestamp                                            |
| ZIP     | [yamaguchi-rh-2026-src.zip](./yamaguchi-rh-2026-src.zip)     | -                                                            |
| ZIP.ots | [yamaguchi-rh-2026-src.zip.ots](./yamaguchi-rh-2026-src.zip.ots) | Bitcoin timestamp                                            |
| Source  | [src/](./src/)                                               | 77 .c programs, 15 .h headers                                |

---

## License

| Content             | License          |
| ------------------- | ---------------- |
| Source code (.c/.h) | **AGPLv3**       |
| Article (.tex/.pdf) | **CC-BY-SA 4.0** |

---

## Result

The paper constructs the **Gram Jacobi matrix** $J_N$ — a finite-dimensional
Hermitian matrix whose eigenvalues approximate the imaginary parts $\gamma_n$ of
the non-trivial zeros of $\zeta(s)$ — and proves the **Riemann Hypothesis**.

**Proof chain (§5):**

- **Lemma I:** Gram Jacobi construction — self-adjoint, Weyl law density $\rho(E) = \frac{1}{2\pi}\log\frac{E}{2\pi}$ (§6, Theorem 6.1, Appendix C)
- **Lemma II:** Correction formula $\delta a_n = -\pi(S(\gamma_n^+) - 0.5)/\theta'(g_{n-1})$, RMS 0.0090 (§6.3, Lemma 6.2, Appendices B–C)
- **Lemma III:** Sturm oscillation $\arg\det(J_N - EI) = \pi N_J(E)$ (§6, Lemma 6.3, Appendix D)
- **Convergence:** Abel–Fejér–Moore–Osgood + Guinand–Weil explicit formula (§9)
- **Theorem I:** Trace formula $\operatorname{Tr} h(J_\infty) = \sum_k h(\gamma_k)$ for all Schwartz $h$, via Birman–Krein + Riemann–von Mangoldt IBP + Guinand–Weil (§10.3)
- **Theorem II:** Spectral determinant identity $D_N(z)/\xi(\tfrac12+iz) \to c$ via $\det_2$ Hilbert–Schmidt regularisation (§10.1, Appendix E)
- **Theorem III:** Hadamard rigidity — evenness + ratio convergence $\implies \{\pm\lambda_k\} = \{\pm\gamma_k\}$ (§10.1)
- **RH:** Self-adjointness of $J_\infty$ $\implies$ all eigenvalues real $\implies$ all zeros on $\Re(s) = \tfrac12$ (§10.4)

---

## Key Numerical Results

| Metric | Value |
|--------|-------|
| Correction formula RMS | 0.0090 (99.9% variance explained) |
| Heat kernel trace ratio | 0.9999996 |
| O(1/√N) error bound | RMS∞ = 0.61 |
| Killip–Simon sum rule | Σ(b/a)² = 0.059 < 1 |
| Geronimo–Case smoothness | Σ\|Δ(b/a)\| = 0.069 < ∞ |
| Carleman condition | Σ 1/b_k ~ 0.4√N → ∞ |
| Hadamard product | verified at 10,000 zeros (mpmath, 18-digit) |
| D_N/ξ ratio | → c ≈ 1.96 (Hadamard rigidity) |
| Prime detection | 25/25 primes p ≤ 97 via spectral shift DFT |
| Contradiction machine | 80/80 falsified (SNR > 2.6 × 10⁹) |
| de Boor–Golub forward check | < 10⁻¹³ |

---

## Source Code (src/)

77 proof-essential programs, 15 reference data headers. Build: `make`

### Lemma I — Gram Jacobi construction (§6)

| File | Purpose |
|------|---------|
| `deboor_golub.c` | de Boor–Golub inverse reconstruction (forward check < 10⁻¹³) |
| `weyl_law_verify.c` | Weyl law, Killip–Simon (0.059), Geronimo–Case (0.069) |
| `phase_shooting.c` | Prufer phase shooting for eigenvalue shift |
| `stronger_conditions.c` | Carleman condition, Janas–Naboko, Schatten class |
| `isospectral_flow.c` | Isospectral manifold flow |
| `trace_class.c` | Killip–Simon trace-class verification |
| `jacobian_analysis.c` | Jacobian of eigenvalue map |
| `rescaled_operator.c` | Rescaled operator N → ∞ |
| `gauge_invariance.c` | Midpoint gauge as N → ∞ limit |
| `block_jacobi.c` | Block Jacobi construction |
| `analytic_detrend.c` | Detrended spectral density |
| `analytic_entry_formula.c` | Analytic entry formula verification |
| `det_xi_match.c` | Eigenvalue counting function vs RvM comparison |
| `random_matrix_ou.c` | Random matrix OU process: GOE to Jacobi |
| `forward_prime.c` | Forward construction: primes → Jacobi → eigenvalues |
| `explicit_test.c` | Explicit formula test: spectral shift vs prime sum |

### Lemma II — Correction formula (§6.3)

| File | Purpose |
|------|---------|
| `derive_k.c` | Correction formula verification (RMS 0.0090, Pearson 0.9997) |
| `derive_k2.c` | Second-order analysis (linear optimal, +0.16% only) |
| `derive_k_gmp.c` | 333-bit GMP precision verification |

### Convergence — Abel–Fejér–Moore–Osgood (§9)

| File | Purpose |
|------|---------|
| `classify_primes.c` | Locked vs outlier prime classification |
| `kp_verify.c` | k(p) = A + B log(p) verification |
| `derive_kp.c` | Analytic k(p) derivation |
| `predict.c` | Prime-power hierarchy A_m = A₁^m/m (p ≤ 47) |
| `sincos_n25.c` | Fejér DFT, sin/cos channel decomposition |
| `spectral_shift.c` | Spectral shift DFT (detects 25/25 primes) |
| `prime_detector.c` | 100% target matrix as prime-frequency detector |
| `prime_resonance.c` | Map primes to eigenvalues via Feynman–Hellmann |
| `harmonic_test.c` | Test if the 20% gap is prime-power harmonics |
| `frac.c` | Fractional prime-power harmonic search |
| `fractional_test.c` | Fractional frequency search at half-integer exponents |
| `legendre_test.c` | Legendre symbol matrices for prime-correlated eigenvalues |
| `test_epsilon_paths.c` | Epsilon → 0 path independence of DFT limit |
| `test_fejer_prime_sum.c` | Fejér kernel convergence for p⁻¹/² divergence |

### Theorem I — Trace formula (§10.3)

| File | Purpose |
|------|---------|
| `trace_verify.c` | Heat kernel trace + moment traces (ratio 0.9999996) |
| `trace_error_bound.c` | O(1/√N) error bound verification |
| `heat_kernel_expansion.c` | Heat kernel asymptotics for Gram Jacobi |
| `tauberian_argument.c` | Uniform Tauberian theorem verification |
| `stationary_phase_test.c` | Kusmin–Landau stationary-phase bound |
| `phase234.c` | Phase 2–4: resolvent, m-function, convergence |
| `potential_search.c` | Custom potential for RvM density matching |
| `lanczos_potential.c` | Lanczos-accelerated large-N Dyson BM |

### Theorem II — Spectral determinant (§10.1, §11)

| File | Purpose |
|------|---------|
| `det2_vs_xi.c` | D_N/ξ ratio convergence → c ≈ 1.96 |
| `det_uniform_bound.c` | Determinant uniform bound |
| `det2_uniform_bound.c` | D_N uniform bound |
| `spectral_rigidity.c` | D_N evenness verification (machine precision) |
| `hadamard_vs_analytic.c` | Hadamard product vs analytic xi(s) |
| `hadamard_extrapolation.c` | Richardson extrapolation for Hadamard tail |
| `hadamard_10k.c` | Hadamard product at 10,000 zeros (18-digit) |
| `hadamard_terminal.c` | Contradiction machine (80/80, SNR > 2.6 × 10⁹) |
| `xi_hadamard_vs_mpmath.c` | 18-digit mpmath comparison |

### Theorem III — Hadamard rigidity (§10.1)

| File | Purpose |
|------|---------|
| `spectral_chain.c` | Unified spectral proof chain |
| `spectral_chain_extended.c` | Spectral chain extended (Euler form) |

### Paths A/B/C — Independent verification (§11.8)

| File | Purpose |
|------|---------|
| `prove_path_a_determinant.c` | Birman–Krein spectral determinant path |
| `prove_path_b_gaussian.c` | Gaussian-regularised Guinand–Weil path |
| `prove_epsilon_zero_closure.c` | Fejér kernel closure proof: ε → 0 |

### Closures and m-function bridge

| File | Purpose |
|------|---------|
| `l8_verify.c` | WKB resolvent phase C₀ = 1.3286 |
| `resolvent_trace.c` | Resolvent trace convergence |
| `krein_ssf_enhanced.c` | Enhanced Krein SSF verification |
| `mfunction_bridge.c` | m-function bridge: m(z) vs zeta'/zeta |
| `m_function.c` | Weyl m-function via continued fraction |
| `m_zeta_bridge.c` | m–ζ bridge (dispersion relation) |
| `spectral_cauchy.c` | Spectral Cauchy transform |
| `unnormalized_cauchy.c` | Unnormalized Cauchy transform |
| `riemann_siegel_zeta.c` | Riemann–Siegel Z(t) via Euler–Maclaurin |
| `wkb_resolvent.c` | WKB resolvent phase analysis |

### Isospectral gauge exploration (§7)

| File | Purpose |
|------|---------|
| `prime_fit.c` | Optimize isospectral gauge for prime perturbation |
| `optimal_jacobi.c` | Reconstruct optimal Jacobi from least-squares gauge |
| `complex_sym.c` | Complex symmetric Jacobi via gauge potential |
| `dual_gauge.c` | Dual independent gauge spectra for sin + cos channels |
| `full_gauge_100.c` | All gauge mechanisms combined (21 params) |
| `gauge_fit_100.c` | Gauge potential sufficiency test |
| `final_optimum.c` | Combined optimization on all 25 gauge parameters |
| `explicit_optimum.c` | Coordinate descent on 25D gauge manifold |
| `block_dbg.c` | 2×2 block dBG with diagonal weight matrices |
| `block_deboor.c` | 2×2 block dBG with diagonal coupling |
| `block_coupled.c` | True coupled 2×2 block dBG (Younes–Trogdon) |
| `block_fit_100.c` | Dual-channel dBG for 100% R² |

### Reference data headers

`refdata_100.h`, `refdata_200.h`, `refdata_1000.h`, `refdata_2000.h`, `zetazero_refdata10000.h`, `random_matrix_utils.h`, `arb.h`

### Specialized reference data

`refdata_c13.h`, `refdata_c13_N5_T200.h`, `refdata_c13_N20_T200.h`, `refdata_c13_N20_T200_dps100.h`, `refdata_c13_N30_T200.h`, `refdata_c17_N5_T200.h`, `refdata_c30_N5_T200.h`, `refdata_c50_N5_T200.h`
