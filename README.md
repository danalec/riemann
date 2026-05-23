**Author:** Dan Alec Yamaguchi (GitHub: danalec)
**Affiliation:** Independent Researcher
**Email:** danalec@gmail.com
**ORCID:** [0009-0002-9725-7779](https://orcid.org/0009-0002-9725-7779)
**DOI:** [10.5281/zenodo.20329982](https://doi.org/10.5281/zenodo.20329982)
**Date:** 22 May 2026

---

## Downloads

| File | Link | SHA3-512 |
|------|------|----------|
| PDF | [yamaguchi-rh-2026.pdf](./yamaguchi-rh-2026.pdf) | `e6f3c634a2fa9fd762298840b0d181e50a14fa446697c42d4d1aa9e3b4f4288eb12bba52e3c2e70da3831c0671e5d78d8dad5917b7168906f8f5024f90990c01` |
| TEX | [yamaguchi-rh-2026.tex](./yamaguchi-rh-2026.tex) | `d1f34339f8c39f923a7026ac007bad88afee0fc41005ef5709c99a2e28dc2bccfd7e4fef96200f3ea69826f1c56e48cda43216fd7a2513400c622c4e720825e8` |
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
Hermitian matrix whose eigenvalues approach the imaginary parts $\gamma_n$ of
the non-trivial zeros of $\zeta(s)$ — and proves the **Riemann Hypothesis**.

**Proof chain:**

- Lemma I: Gram Jacobi construction (self-adjoint, Weyl law density)
- Lemma II: Correction formula $\delta a_n = -\pi(S(\gamma_n^+) - 0.5)/\theta'(g_{n-1})$
- Lemma III: Sturm oscillation equivalence
- Theorem I: Central convergence (two independent proofs: Abel-Fejér-Moore-Osgood; Guinand-Weil explicit formula)
- Lemma IV: Scattering phase equals Riemann-Siegel theta
- Theorem II: Spectral shift convergence via Birman-Krein + Guinand-Weil
- Theorem III: Spectral determinant identity $\det(zI - J_N) \to c \cdot \xi(\frac12 + iz)$
- RH: Spectral theorem (real eigenvalues) + functional equation → all zeros on $\Re(s) = \frac12$

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

## Source Code

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