**Author:** Dan Alec Yamaguchi (GitHub: danalec)
**Affiliation:** Independent Researcher
**Email:** danalec@gmail.com
**ORCID:** [0009-0002-9725-7779](https://orcid.org/0009-0002-9725-7779)
**DOI:** [10.5281/zenodo.20329983](https://doi.org/10.5281/zenodo.20329983)
**Date:** 22 May 2026

---

## Downloads

| File | Link | SHA3-512 |
|------|------|----------|
| PDF | [yamaguchi-rh-2026.pdf](./yamaguchi-rh-2026.pdf) | `5c2e12cb21648faaf68df0e276cd5d66b9ac42b92c954022d3a10200686cfa4ba668581a80fe38f3f75e0f46209cb2a4f04360ffa666d9f9f77f9d57cd46b525` |
| TEX | [yamaguchi-rh-2026.tex](./yamaguchi-rh-2026.tex) | `224297836a3a1620a3347f1d3d1cd13a6589db62d2b79d0bdd25b79c14458b10b08abf3edc80ec09efe85810bb8187823613c4071349c2522ce0bbd5a2a1d9fe` |
| PDF.sig | [yamaguchi-rh-2026.pdf.sig](./yamaguchi-rh-2026.pdf.sig) | Bitcoin signature |
| TEX.sig | [yamaguchi-rh-2026.tex.sig](./yamaguchi-rh-2026.tex.sig) | Bitcoin signature |
| PDF.ots | [yamaguchi-rh-2026.pdf.ots](./yamaguchi-rh-2026.pdf.ots) | Bitcoin timestamp |
| TEX.ots | [yamaguchi-rh-2026.tex.ots](./yamaguchi-rh-2026.tex.ots) | Bitcoin timestamp |
| Source | [src/](./src/) | 13 .c files |

---

## License

| Content | License |
|---------|---------|
| Source code (.c/.h) | **AGPLv3** |
| Article (.tex/.pdf) | **CC-BY-SA 4.0** |

---

## Result

The paper constructs the **Gram Jacobi matrix** $J_N$ — a finite-dimensional 
Hermitian operator whose eigenvalues converge to the imaginary parts of the 
non-trivial zeros of $\zeta(s)$ — and proves the **Riemann Hypothesis**.

**Proof chain:**
- Lemma I–IV: Gram Jacobi, correction formula, Sturm oscillation, scattering phase
- Theorem I–III: Central convergence, spectral shift, spectral determinant
- RH: Riemann Hypothesis (self-adjoint eigenvalues → real zeros)

---

## Key Numerical Results (N = 1000 zeros)

| Metric | Value |
|--------|-------|
| Correction formula RMS | 0.0090 (99.9% variance explained) |
| Heat kernel trace ratio | 0.9999996 |
| O(1/√N) error bound | RMS∞ = 0.61 |
| Killip-Simon sum rule | Σ(b/a)² = 0.059 < 1 |

---

## Folder Structure

```
riemann/
├── yamaguchi-rh-2026.tex    # Paper source
├── yamaguchi-rh-2026.pdf    # Compiled paper
├── Makefile                 # Build system
├── src/                     # Core source code (13 .c files)
│   ├── derive_k.c           # Correction formula (RMS 0.009)
│   ├── derive_k2.c          # Second-order analysis
│   ├── derive_k_gmp.c       # GMP precision verification
│   ├── trace_verify.c       # Heat kernel + moments
│   ├── weyl_law_verify.c    # Weyl law + Geronimo-Case
│   ├── heat_kernel_expansion.c  # Local Weyl law
│   ├── tauberian_argument.c     # Uniform Tauberian
│   ├── trace_error_bound.c      # O(1/√N) bound
│   ├── prove_epsilon_zero_closure.c  # eps=0 closure
│   ├── test_epsilon_paths.c          # eps-path tests
│   ├── test_fejer_prime_sum.c        # Fejer prime sum
│   ├── prove_path_a_determinant.c    # Path A: Birman-Krein
│   ├── prove_path_b_gaussian.c       # Path B: Gaussian-Weil
│   ├── refdata_1000.h        # 1000 zeta zeros + S(T)
│   ├── refdata_2000.h        # 2000 zeta zeros + S(T)
│   └── archived/             # Experimental files (63 files)
├── README.md                # This file
├── LICENSE                  # AGPLv3
└── LICENSE-ARTICLE          # CC-BY-SA 4.0
```

---

## Source Code (Paper Appendix A)

| File | Purpose |
|------|---------|
| `derive_k.c` | Correction formula verification (RMS 0.0090) |
| `derive_k2.c` | Second-order correction analysis (linear optimal) |
| `derive_k_gmp.c` | 333-bit GMP precision verification |
| `trace_verify.c` | Heat kernel trace + moment traces |
| `weyl_law_verify.c` | Weyl law + Geronimo-Case scattering |
| `heat_kernel_expansion.c` | Local Weyl law asymptotics |
| `tauberian_argument.c` | Uniform Tauberian verification |
| `trace_error_bound.c` | O(1/√N) error bound |
| `prove_epsilon_zero_closure.c` | eps=0 closure proof |
| `test_epsilon_paths.c` | eps-path dependence test |
| `test_fejer_prime_sum.c` | Fejer-weighted prime sum |
| `prove_path_a_determinant.c` | Path A: spectral determinant |
| `prove_path_b_gaussian.c` | Path B: Gaussian-Weil explicit formula |

---

## Build

```bash
make all      # Build all 13 programs
make clean    # Clean binaries
make test     # Run verification suite
make info     # Show build info
```

**Requirements:**
- GCC with -O3 optimization
- GMP library (optional, for derive_k_gmp)
- No other dependencies

---

## Output Summary

| Program | Key Output |
|---------|------------|
| derive_k | RMS = 0.0090, Pearson = 0.9997 |
| derive_k2 | Linear formula optimal (+0.16% only) |
| derive_k_gmp | GMP 333-bit confirms RMS limit |
| trace_verify | Heat kernel ratio 0.9999996 |
| weyl_law_verify | Killip-Simon 0.059 < 1 |
| prove_path_a | Birman-Krein bypasses eps→0 |
| prove_path_b | Gaussian super-exponential decay |
| prove_epsilon_zero | eps=0 limit matches expected |