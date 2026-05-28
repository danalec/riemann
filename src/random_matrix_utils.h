/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * random_matrix_utils.h - Shared utilities
 *
 * Provides:
 *   - Complex arithmetic via Euler exponential form:
 *     sin(z) = (e^{iz} - e^{-iz})/2i,  cos(z) = (e^{iz} + e^{-iz})/2
 *     Every prime-frequency oscillation is expressed as p^{ik} = e^{ik ln p}
 *   - xmalloc, PRNG seeding, Jacobi eigensolver
 */

#pragma once

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.71828182845904523536
#endif

__attribute__((unused)) static void *xmalloc(size_t sz) {
  void *p = malloc(sz);
  if (!p) {
    fprintf(stderr, "xmalloc(%llu) failed\n", (unsigned long long)sz);
    exit(1);
  }
  return p;
}

static inline void *xcalloc(size_t nmemb, size_t sz) {
  void *p = calloc(nmemb, sz);
  if (!p) {
    fprintf(stderr, "xcalloc(%llu,%llu) failed\n", (unsigned long long)nmemb,
            (unsigned long long)sz);
    exit(1);
  }
  return p;
}

/* ==========================================================================
 *   xoshiro256** PRNG — full
 * ========================================================================== */

static uint64_t rng_s[4];

__attribute__((unused)) static uint64_t splitmix64(uint64_t *state) {
  uint64_t z = (*state += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

__attribute__((unused)) static void rng_set_seed(uint64_t seed) {
  uint64_t st = seed;
  rng_s[0] = splitmix64(&st);
  rng_s[1] = splitmix64(&st);
  rng_s[2] = splitmix64(&st);
  rng_s[3] = splitmix64(&st);
}

static inline uint64_t rng_next64(void) {
  uint64_t *s = rng_s;
  uint64_t result = s[0] * 0x9e3779b97f4a7c15ULL;
  uint64_t t = s[1] << 17;
  s[2] ^= s[0];
  s[3] ^= s[1];
  s[1] ^= s[2];
  s[0] ^= s[3];
  s[2] ^= t;
  s[3] = (s[3] << 45) | (s[3] >> 19);
  return result;
}

static inline double rng_uniform(void) {
  return (double)(rng_next64() >> 11) * 0x1.0p-53;
}

static inline double rng_normal(void) {
  double u1 = rng_uniform();
  double u2 = rng_uniform();
  while (u1 < 1e-15)
    u1 = rng_uniform();
  return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/* ==========================================================================
 *   Matrix utilities (row-major, N×N)
 * ========================================================================== */

__attribute__((unused)) static inline void mat_copy(double *dst,
                                                    const double *src, int N) {
  memcpy(dst, src, (size_t)(N * N) * sizeof(double));
}

__attribute__((unused)) static inline double mat_trace(const double *M, int N) {
  double tr = 0.0;
  for (int i = 0; i < N; i++)
    tr += M[i * N + i];
  return tr;
}

__attribute__((unused)) static inline double mat_frobenius_sq(const double *M,
                                                              int N) {
  double s = 0.0;
  for (int i = 0; i < N * N; i++)
    s += M[i] * M[i];
  return s;
}

__attribute__((unused)) static inline void mat_add_goe_noise(double *M, int N,
                                                             double sigma) {
  double s2 = sigma / sqrt(2.0 * (double)N);
  for (int i = 0; i < N; i++) {
    M[i * N + i] += s2 * rng_normal();
    for (int j = i + 1; j < N; j++) {
      double g = s2 * rng_normal();
      M[i * N + j] += g;
      M[j * N + i] += g;
    }
  }
}

/* ==========================================================================
 *   Jacobi eigenvalue decomposition (no eigenvectors)
 *   Classical cyclic Jacobi for real symmetric matrices.
 * ========================================================================== */

__attribute__((unused)) static inline void
jacobi_eigenvalues(const double *A, int N, double *evals) {
  double *V = (double *)xmalloc((size_t)(N * N) * sizeof(double));
  memcpy(V, A, (size_t)(N * N) * sizeof(double));

  for (int sweep = 0; sweep < 50; sweep++) {
    double max_off = 0.0;
    for (int p = 0; p < N - 1; p++)
      for (int q = p + 1; q < N; q++) {
        double v = fabs(V[p * N + q]);
        if (v > max_off)
          max_off = v;
      }
    if (max_off < 1e-14)
      break;

    for (int p = 0; p < N - 1; p++) {
      for (int q = p + 1; q < N; q++) {
        double apq = V[p * N + q];
        if (fabs(apq) < 1e-16 * (fabs(V[p * N + p]) + fabs(V[q * N + q]) + 1.0))
          continue;
        double app = V[p * N + p], aqq = V[q * N + q];
        double tau = (aqq - app) / (2.0 * apq);
        double t;
        if (tau >= 0.0)
          t = 1.0 / (tau + sqrt(1.0 + tau * tau));
        else
          t = -1.0 / (-tau + sqrt(1.0 + tau * tau));
        double c = 1.0 / sqrt(1.0 + t * t);
        double s = t * c;
        for (int i = 0; i < N; i++) {
          double vip = V[i * N + p], viq = V[i * N + q];
          V[i * N + p] = vip * c - viq * s;
          V[i * N + q] = vip * s + viq * c;
        }
        for (int j = 0; j < N; j++) {
          double vpj = V[p * N + j], vqj = V[q * N + j];
          V[p * N + j] = vpj * c - vqj * s;
          V[q * N + j] = vpj * s + vqj * c;
        }
      }
    }
  }

  for (int i = 0; i < N; i++)
    evals[i] = V[i * N + i];
  for (int i = 1; i < N; i++) {
    double key = evals[i];
    int j = i - 1;
    while (j >= 0 && evals[j] > key) {
      evals[j + 1] = evals[j];
      j--;
    }
    evals[j + 1] = key;
  }
  free(V);
}

/* ==========================================================================
 *   Euler Complex Arithmetic
 *
 *   Fundamental relation:  sin(z) = (e^{iz} - e^{-iz}) / 2i
 *                          cos(z) = (e^{iz} + e^{-iz}) / 2
 *
 *   For prime-frequency analysis, every oscillation is a complex exponential:
 *     sin(k ln p) = [p^{ik} - p^{-ik}] / 2i
 *     cos(k ln p) = [p^{ik} + p^{-ik}] / 2
 *
 *   This makes the TWO-SIDED spectrum explicit:
 *   each prime p contributes conjugate pair at frequencies ±ln(p).
 *
 *   The prime perturbation α_p sin(k ln p) becomes:
 *     α_p/(2i) · p^{ik}  +  (-α_p/(2i)) · p^{-ik}
 *   with amplitudes  c_p⁺ = iα_p/2  and  c_p⁻ = -iα_p/2 = conj(c_p⁺)
 * ========================================================================== */

typedef struct {
  double re, im;
} Cpx;

static inline Cpx cpx_make(double re, double im) {
  Cpx r;
  r.re = re;
  r.im = im;
  return r;
}
static inline Cpx cpx_add(Cpx a, Cpx b) {
  return cpx_make(a.re + b.re, a.im + b.im);
}
static inline Cpx cpx_sub(Cpx a, Cpx b) {
  return cpx_make(a.re - b.re, a.im - b.im);
}
static inline Cpx cpx_mul(Cpx a, Cpx b) {
  return cpx_make(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}
static inline Cpx cpx_scale(double s, Cpx a) {
  return cpx_make(s * a.re, s * a.im);
}
static inline double cpx_abs(Cpx a) { return sqrt(a.re * a.re + a.im * a.im); }
static inline double cpx_abs2(Cpx a) { return a.re * a.re + a.im * a.im; }
static inline Cpx cpx_conj(Cpx a) { return cpx_make(a.re, -a.im); }
static inline Cpx cpx_neg(Cpx a) { return cpx_make(-a.re, -a.im); }
static inline Cpx cpx_inv(Cpx a) {
  double d = cpx_abs2(a);
  return cpx_make(a.re / d, -a.im / d);
}
static inline Cpx cpx_div(Cpx a, Cpx b) { return cpx_mul(a, cpx_inv(b)); }
static inline double cpx_arg(Cpx a) { return atan2(a.im, a.re); }
static inline Cpx cpx_sqrt(Cpx a) {
  double r = cpx_abs(a);
  double t = cpx_arg(a);
  return cpx_make(sqrt(r) * cos(t / 2), sqrt(r) * sin(t / 2));
}

/* e^{ix} = cos(x) + i sin(x) — the fundamental building block */
static inline Cpx cpx_exp_i(double x) { return cpx_make(cos(x), sin(x)); }

/* e^{-ix} = cos(x) - i sin(x) */
static inline Cpx cpx_exp_mi(double x) { return cpx_make(cos(x), -sin(x)); }

/* p^{ik} = e^{ik ln p} — prime-frequency complex exponential */
static inline Cpx cpx_p_ik(double p, double k) { return cpx_exp_i(k * log(p)); }

/* p^{-ik} = e^{-ik ln p} — conjugate frequency */
static inline Cpx cpx_p_mik(double p, double k) {
  return cpx_exp_mi(k * log(p));
}

/* Euler form: sin(k ln p) = [p^{ik} - p^{-ik}] / 2i */
static inline double euler_sin(double p, double k) {
  Cpx pos = cpx_p_ik(p, k), neg = cpx_p_mik(p, k);
  return (pos.im - neg.im) / 2.0;
}

/* Euler form: cos(k ln p) = [p^{ik} + p^{-ik}] / 2 */
static inline double euler_cos(double p, double k) {
  Cpx pos = cpx_p_ik(p, k), neg = cpx_p_mik(p, k);
  return (pos.re + neg.re) / 2.0;
}

/* Euler form: complex amplitude pair for α_p sin(k ln p)
 * Returns c_p⁺ = iα_p/2  (coefficient of p^{ik})
 * Its conjugate conj(c_p⁺) = -iα_p/2 is the coefficient of p^{-ik} */
static inline Cpx euler_alpha_coeff(double p) {
  double alpha = -log(p) / (2.0 * M_PI * sqrt(p));
  return cpx_make(0.0, alpha / 2.0);
}

/* Full prime perturbation at site k via Euler decomposition:
 * δ(k) = Re[ Σ_p 2·c_p⁺ · p^{ik} ]
 *       = Σ_p [ Re(c_p⁺)·cos(k ln p) - Im(c_p⁺)·sin(k ln p) ] × 2
 *       = Σ_p α_p · sin(k ln p)   [when c_p⁺ = iα_p/2] */
static inline double euler_prime_sum(const int *primes, int np, double k) {
  Cpx sum = cpx_make(0, 0);
  for (int pi = 0; pi < np; pi++) {
    double p = (double)primes[pi];
    Cpx cp = euler_alpha_coeff(p);
    Cpx phase = cpx_p_ik(p, k);
    sum = cpx_add(sum, cpx_mul(cp, phase));
  }
  return 2.0 * sum.re;
}

/* DFT at a specific frequency ω:  (1/N) Σ x_k · e^{-iωk}
 * This directly extracts the complex amplitude of the e^{iωk} component */
static inline Cpx dft_at_freq(const double *x, int N, double omega) {
  Cpx sum = cpx_make(0, 0);
  for (int k = 0; k < N; k++) {
    Cpx twiddle = cpx_exp_i(-omega * (double)k);
    sum = cpx_add(sum, cpx_scale(x[k], twiddle));
  }
  sum.re /= (double)N;
  sum.im /= (double)N;
  return sum;
}

/* DFT of complex array at frequency ω */
static inline Cpx dft_complex_at_freq(const Cpx *x, int N, double omega) {
  Cpx sum = cpx_make(0, 0);
  for (int k = 0; k < N; k++) {
    Cpx twiddle = cpx_exp_i(-omega * (double)k);
    sum = cpx_add(sum, cpx_mul(x[k], twiddle));
  }
  sum.re /= (double)N;
  sum.im /= (double)N;
  return sum;
}

/* Conjugate symmetry check: DFT[-ω] should = conj(DFT[+ω]) for real signals */
static inline double conj_symmetry_error(const double *x, int N, double omega) {
  Cpx dp = dft_at_freq(x, N, omega);
  Cpx dm = dft_at_freq(x, N, -omega);
  Cpx expected_conj = cpx_conj(dp);
  return cpx_abs(cpx_sub(dm, expected_conj));
}
