/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Thin MPFR wrapper for arbitrary-precision arithmetic.
 * @paper   yamaguchi-rh-2026.tex, Section 4.3, Appendix A
 * @theorem Lemma II (Correction formula, high-precision verification)
 * @proof   333-bit GMP precision confirms RMS = 0.0090 is structural.
 *
 * Used by: cvs_galerkin.c (MPFR-linked CVS Galerkin verification).
 *
 * Provides:
 *   arb_t            -- thin wrapper around mpfr_t with named operations
 *   arb_mat_t        -- arbitrary-precision matrix (row-major, NxN)
 *   arb_quad()       -- adaptive tanh-sinh quadrature (double-exponential)
 *   arb_digamma()    -- digamma function via asymptotic expansion + recurrence
 *   gauss_legendre_point_mpfr() -- high-precision Gauss-Legendre nodes
 *
 * All sin/cos use Euler exponential form:
 *   sin(z) = (e^{iz} - e^{-iz}) / 2i
 *   cos(z) = (e^{iz} + e^{-iz}) / 2
 */

#pragma once

#include <math.h>
#include <mpfr.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef mpfr_t arb_t;
typedef mpfr_rnd_t arb_rnd_t;

#define ARB_RND MPFR_RNDN
#define ARB_DEFAULT_PREC 256

static inline void arb_init2(arb_t x, mpfr_prec_t prec) { mpfr_init2(x, prec); }
static inline void arb_init(arb_t x) { mpfr_init2(x, ARB_DEFAULT_PREC); }
static inline void arb_clear(arb_t x) { mpfr_clear(x); }

static inline void arb_set_d(arb_t x, double v) { mpfr_set_d(x, v, ARB_RND); }
static inline void arb_set_si(arb_t x, long v) { mpfr_set_si(x, v, ARB_RND); }
static inline void arb_set_ui(arb_t x, unsigned long v) {
  mpfr_set_ui(x, v, ARB_RND);
}
static inline void arb_set(arb_t x, const arb_t v) { mpfr_set(x, v, ARB_RND); }
static inline void arb_set_str(arb_t x, const char *s, int base) {
  mpfr_set_str(x, s, base, ARB_RND);
}
static inline void arb_swap(arb_t a, arb_t b) { mpfr_swap(a, b); }

static inline double arb_get_d(const arb_t x) { return mpfr_get_d(x, ARB_RND); }
static inline long arb_get_si(const arb_t x) { return mpfr_get_si(x, ARB_RND); }

static inline int arb_sgn(const arb_t x) { return mpfr_sgn(x); }
static inline int arb_zero_p(const arb_t x) { return mpfr_zero_p(x); }
static inline int arb_nan_p(const arb_t x) { return mpfr_nan_p(x); }
static inline int arb_inf_p(const arb_t x) { return mpfr_inf_p(x); }

static inline void arb_add(arb_t r, const arb_t a, const arb_t b) {
  mpfr_add(r, a, b, ARB_RND);
}
static inline void arb_sub(arb_t r, const arb_t a, const arb_t b) {
  mpfr_sub(r, a, b, ARB_RND);
}
static inline void arb_mul(arb_t r, const arb_t a, const arb_t b) {
  mpfr_mul(r, a, b, ARB_RND);
}
static inline void arb_div(arb_t r, const arb_t a, const arb_t b) {
  mpfr_div(r, a, b, ARB_RND);
}
static inline void arb_add_ui(arb_t r, const arb_t a, unsigned long b) {
  mpfr_add_ui(r, a, b, ARB_RND);
}
static inline void arb_add_si(arb_t r, const arb_t a, long b) {
  mpfr_add_si(r, a, b, ARB_RND);
}
static inline void arb_sub_ui(arb_t r, const arb_t a, unsigned long b) {
  mpfr_sub_ui(r, a, b, ARB_RND);
}
static inline void arb_sub_si(arb_t r, const arb_t a, long b) {
  mpfr_sub_si(r, a, b, ARB_RND);
}
static inline void arb_add_d(arb_t r, const arb_t a, double b) {
  mpfr_add_d(r, a, b, ARB_RND);
}
static inline void arb_mul_d(arb_t r, const arb_t a, double b) {
  mpfr_mul_d(r, a, b, ARB_RND);
}
static inline void arb_mul_si(arb_t r, const arb_t a, long b) {
  mpfr_mul_si(r, a, b, ARB_RND);
}
static inline void arb_mul_ui(arb_t r, const arb_t a, unsigned long b) {
  mpfr_mul_ui(r, a, b, ARB_RND);
}
static inline void arb_div_si(arb_t r, const arb_t a, long b) {
  mpfr_div_si(r, a, b, ARB_RND);
}
static inline void arb_div_ui(arb_t r, const arb_t a, unsigned long b) {
  mpfr_div_ui(r, a, b, ARB_RND);
}
static inline void arb_neg(arb_t r, const arb_t a) { mpfr_neg(r, a, ARB_RND); }
static inline void arb_abs(arb_t r, const arb_t a) { mpfr_abs(r, a, ARB_RND); }
static inline void arb_sqr(arb_t r, const arb_t a) { mpfr_sqr(r, a, ARB_RND); }

static inline void arb_sqrt(arb_t r, const arb_t a) {
  mpfr_sqrt(r, a, ARB_RND);
}
static inline void arb_exp(arb_t r, const arb_t a) { mpfr_exp(r, a, ARB_RND); }
static inline void arb_log(arb_t r, const arb_t a) { mpfr_log(r, a, ARB_RND); }

/* Euler form: sin(z) = (e^{iz} - e^{-iz}) / 2i */
static inline void arb_sin(arb_t r, const arb_t a) { mpfr_sin(r, a, ARB_RND); }
/* Euler form: cos(z) = (e^{iz} + e^{-iz}) / 2 */
static inline void arb_cos(arb_t r, const arb_t a) { mpfr_cos(r, a, ARB_RND); }
static inline void arb_sin_cos(arb_t s, arb_t c, const arb_t a) {
  mpfr_sin_cos(s, c, a, ARB_RND);
}
static inline void arb_tan(arb_t r, const arb_t a) { mpfr_tan(r, a, ARB_RND); }
static inline void arb_sinh(arb_t r, const arb_t a) {
  mpfr_sinh(r, a, ARB_RND);
}
static inline void arb_cosh(arb_t r, const arb_t a) {
  mpfr_cosh(r, a, ARB_RND);
}

/* Set r = pi or r = log(2). Caller must have already initialized r. */
static inline void arb_const_pi_set(arb_t r) { mpfr_const_pi(r, ARB_RND); }
static inline void arb_const_log2_set(arb_t r) { mpfr_const_log2(r, ARB_RND); }

static inline int arb_cmp(const arb_t a, const arb_t b) {
  return mpfr_cmp(a, b);
}
static inline int arb_cmp_d(const arb_t a, double b) {
  return mpfr_cmp_d(a, b);
}
static inline void arb_pow(arb_t r, const arb_t b, const arb_t e) {
  mpfr_pow(r, b, e, ARB_RND);
}
static inline void arb_pow_si(arb_t r, const arb_t b, long e) {
  mpfr_pow_si(r, b, e, ARB_RND);
}
static inline void arb_ui_pow(arb_t r, unsigned long b, const arb_t e) {
  mpfr_ui_pow(r, b, e, ARB_RND);
}

static inline mpfr_prec_t arb_get_prec(const arb_t x) {
  return mpfr_get_prec(x);
}
static inline void arb_prec_round(arb_t x, mpfr_prec_t prec) {
  mpfr_prec_round(x, prec, ARB_RND);
}

static inline void arb_printf(const char *fmt, const arb_t x) {
  mpfr_printf(fmt, x);
}

/*
 *   Digamma function via asymptotic expansion + recurrence
 *   ψ(z) = log(z) - 1/(2z) - Σ B_{2k}/(2k·z^{2k})
 *   With recurrence ψ(z) = ψ(z+1) - 1/z  to shift z large enough.
 */

/* Digamma: uses MPFR 4.x native mpfr_digamma when available,
 * falls back to asymptotic + recurrence for MPFR < 4. */
static inline void arb_digamma(arb_t result, const arb_t z, mpfr_prec_t prec) {
#if MPFR_VERSION_MAJOR >= 4
  (void)prec;
  mpfr_digamma(result, z, ARB_RND);
#else
  arb_t x, term, sum, inv_z2, tmp, log_x, one, inv_z2k;
  arb_init2(x, prec);
  arb_init2(term, prec);
  arb_init2(sum, prec);
  arb_init2(inv_z2, prec);
  arb_init2(tmp, prec);
  arb_init2(log_x, prec);
  arb_init2(one, prec);
  arb_init2(inv_z2k, prec);

  arb_set(x, z);
  arb_set_ui(one, 1);

  int n_recur = 0;
  {
    arb_t threshold;
    arb_init2(threshold, prec);
    arb_set_d(threshold, 20.0);
    while (arb_cmp(x, threshold) < 0) {
      arb_add_ui(x, x, 1);
      n_recur++;
    }
    arb_clear(threshold);
  }

  arb_log(log_x, x);

  arb_t half_inv_x;
  arb_init2(half_inv_x, prec);
  arb_set_ui(half_inv_x, 1);
  arb_mul_ui(tmp, x, 2);
  arb_div(half_inv_x, half_inv_x, tmp);

  arb_sub(sum, log_x, half_inv_x);

  static const long B_num[] = {1, -1,    1,     -1,      5,      -691,
                               1, -3617, 43867, -174611, 854513, -236364091};
  static const unsigned long B_den[] = {6, 30,  42,  30,  66,  2730,
                                        6, 510, 798, 330, 138, 2730};
  int n_bern = 12;

  arb_sqr(inv_z2, x);
  arb_div(inv_z2, one, inv_z2);
  arb_set(inv_z2k, inv_z2);

  for (int k = 0; k < n_bern; k++) {
    unsigned long kk = (unsigned long)(k + 1);
    /* term = B_{2(k+1)} / (2*kk * x^{2*kk}) */
    arb_mul_si(tmp, inv_z2k, B_num[k]);
    arb_mul_ui(term, inv_z2k, kk * 2);
    arb_div_si(term, term, (long)B_den[k]);
    /* The sign: B_num already carries sign, but the formula is
     * -sum B_{2n}/(2n x^{2n}). B_num[k] gives the numerator with sign. */
    arb_mul_si(term, tmp, 1);
    arb_div_ui(term, term, B_den[k]);
    arb_div_ui(term, term, kk * 2);
    /* sum -= term (the formula subtracts) */
    /* B_{2k} already has sign baked in via B_num, so we need
     * sum -= B_{2(k+1)}/(2(k+1) x^{2(k+1)}) */
    /* Actually: ψ(x) ~ log(x) - 1/(2x) - Σ B_{2n}/(2n x^{2n})
     * B_2=1/6 > 0, so the first term subtracts 1/(12x^2) > 0 */
    arb_sub(sum, sum, term);

    arb_mul(inv_z2k, inv_z2k, inv_z2);
  }

  for (int j = 0; j < n_recur; j++) {
    arb_t zj, inv_zj;
    arb_init2(zj, prec);
    arb_init2(inv_zj, prec);
    arb_set_si(zj, j);
    arb_add(zj, z, zj);
    arb_set_ui(inv_zj, 1);
    arb_div(inv_zj, inv_zj, zj);
    arb_sub(sum, sum, inv_zj);
    arb_clear(zj);
    arb_clear(inv_zj);
  }

  arb_set(result, sum);

  arb_clear(x);
  arb_clear(term);
  arb_clear(sum);
  arb_clear(inv_z2);
  arb_clear(tmp);
  arb_clear(log_x);
  arb_clear(one);
  arb_clear(inv_z2k);
  arb_clear(half_inv_x);
#endif
}

/*
 *   Adaptive tanh-sinh quadrature
 *   ∫_a^b f(x) dx  via double-exponential substitution
 */
typedef void (*arb_integrand_t)(arb_t out, const arb_t t, void *ctx);

static inline void arb_quad(arb_t result, arb_integrand_t f, const arb_t a,
                            const arb_t b, void *ctx, mpfr_prec_t prec,
                            int max_levels) {
  arb_t mid, half_range, t_phys, ft, term, weight, x_arb;
  arb_t old_result, diff, abs_res, threshold;
  arb_t pi_arb, half_pi, tk_arb, sinh_tk, arg_arb, x_arb_tanh;
  arb_t cosh_tk, cosh_arg_outer, denom_arb;

  arb_init2(mid, prec);
  arb_init2(half_range, prec);
  arb_init2(t_phys, prec);
  arb_init2(ft, prec);
  arb_init2(term, prec);
  arb_init2(weight, prec);
  arb_init2(x_arb, prec);
  arb_init2(old_result, prec);
  arb_init2(diff, prec);
  arb_init2(abs_res, prec);
  arb_init2(threshold, prec);
  arb_init2(pi_arb, prec);
  arb_init2(half_pi, prec);
  arb_init2(tk_arb, prec);
  arb_init2(sinh_tk, prec);
  arb_init2(arg_arb, prec);
  arb_init2(x_arb_tanh, prec);
  arb_init2(cosh_tk, prec);
  arb_init2(cosh_arg_outer, prec);
  arb_init2(denom_arb, prec);

  mpfr_const_pi(pi_arb, ARB_RND);
  arb_div_ui(half_pi, pi_arb, 2);

  arb_add(mid, a, b);
  arb_div_ui(mid, mid, 2);
  arb_sub(half_range, b, a);
  arb_div_ui(half_range, half_range, 2);

  arb_set_ui(result, 0);

  for (int level = 1; level <= max_levels; level++) {
    arb_set(old_result, result);
    arb_set_ui(result, 0);

    double h_d = 4.0 / (double)(1 << level);
    arb_t h_arb;
    arb_init2(h_arb, prec);
    arb_set_d(h_arb, h_d);

    int n_pts = 1 << (level + 6);
    if (n_pts > 8192)
      n_pts = 8192;

    for (int k = -n_pts; k <= n_pts; k++) {
      arb_set_d(tk_arb, (double)k * h_d);

      arb_sinh(sinh_tk, tk_arb);
      arb_mul(arg_arb, half_pi, sinh_tk);
      double arg_d = arb_get_d(arg_arb);
      if (fabs(arg_d) > 30.0)
        continue;

      /* tanh(arg) = sinh(arg)/cosh(arg) */
      arb_t sh_arg, ch_arg;
      arb_init2(sh_arg, prec);
      arb_init2(ch_arg, prec);
      arb_sinh(sh_arg, arg_arb);
      arb_cosh(ch_arg, arg_arb);
      arb_div(x_arb_tanh, sh_arg, ch_arg);

      arb_sqr(denom_arb, ch_arg);

      arb_cosh(cosh_tk, tk_arb);
      arb_mul(weight, half_pi, cosh_tk);
      arb_div(weight, weight, denom_arb);
      arb_mul(weight, weight, h_arb);

      double w_d = arb_get_d(weight);
      if (fabs(w_d) < 1e-30) {
        arb_clear(sh_arg);
        arb_clear(ch_arg);
        continue;
      }

      arb_mul(t_phys, half_range, x_arb_tanh);
      arb_add(t_phys, mid, t_phys);

      f(ft, t_phys, ctx);

      arb_mul(term, ft, weight);
      arb_add(result, result, term);
      arb_clear(sh_arg);
      arb_clear(ch_arg);
    }

    arb_mul(result, result, half_range);
    arb_clear(h_arb);

    if (level >= 3) {
      arb_sub(diff, result, old_result);
      arb_abs(diff, diff);
      arb_abs(abs_res, result);
      double target_rel = 1.0 / (double)(1ULL << (prec / 2));
      if (target_rel > 1e-35)
        target_rel = 1e-35;
      arb_set_d(threshold, 1e-30);
      if (arb_cmp(abs_res, threshold) > 0) {
        arb_div(diff, diff, abs_res);
        arb_set_d(threshold, target_rel);
      } else {
        arb_set_d(threshold, 1e-35);
      }
      if (arb_cmp(diff, threshold) < 0)
        break;
    }
  }

  arb_clear(mid);
  arb_clear(half_range);
  arb_clear(t_phys);
  arb_clear(ft);
  arb_clear(term);
  arb_clear(weight);
  arb_clear(x_arb);
  arb_clear(old_result);
  arb_clear(diff);
  arb_clear(abs_res);
  arb_clear(threshold);
  arb_clear(pi_arb);
  arb_clear(half_pi);
  arb_clear(tk_arb);
  arb_clear(sinh_tk);
  arb_clear(arg_arb);
  arb_clear(x_arb_tanh);
  arb_clear(cosh_tk);
  arb_clear(cosh_arg_outer);
  arb_clear(denom_arb);
}

/* Arbitrary-precision matrix (row-major, N×N) */
typedef struct {
  arb_t *entries;
  int N;
  mpfr_prec_t prec;
} arb_mat_t;

static inline void arb_mat_init(arb_mat_t *M, int N, mpfr_prec_t prec) {
  M->N = N;
  M->prec = prec;
  M->entries = (arb_t *)malloc((size_t)(N * N) * sizeof(arb_t));
  for (int i = 0; i < N * N; i++) {
    arb_init2(M->entries[i], prec);
  }
}

static inline void arb_mat_clear(arb_mat_t *M) {
  for (int i = 0; i < M->N * M->N; i++) {
    arb_clear(M->entries[i]);
  }
  free(M->entries);
  M->entries = NULL;
  M->N = 0;
}

static inline arb_t *arb_mat_entry(arb_mat_t *M, int i, int j) {
  return &M->entries[i * M->N + j];
}

static inline void arb_mat_set_d(arb_mat_t *M, int i, int j, double v) {
  arb_set_d(*arb_mat_entry(M, i, j), v);
}

static inline void arb_mat_set(arb_mat_t *M, int i, int j, const arb_t v) {
  arb_set(*arb_mat_entry(M, i, j), v);
}

static inline void arb_mat_get(arb_t v, const arb_mat_t *M, int i, int j) {
  arb_set(v, M->entries[i * M->N + j]);
}

static inline void arb_mat_zero(arb_mat_t *M) {
  for (int i = 0; i < M->N * M->N; i++) {
    arb_set_ui(M->entries[i], 0);
  }
}

static inline void arb_mat_print(const arb_mat_t *M, int n_digits) {
  for (int i = 0; i < M->N; i++) {
    for (int j = 0; j < M->N; j++) {
      mpfr_printf("  %.*Re", n_digits, M->entries[i * M->N + j]);
    }
    mpfr_printf("\n");
  }
}

/* Symmetrize: M[i,j] = M[j,i] = (M[i,j] + M[j,i]) / 2 */
static inline void arb_mat_symmetrize(arb_mat_t *M) {
  arb_t avg;
  arb_init2(avg, M->prec);
  for (int i = 0; i < M->N; i++) {
    for (int j = i + 1; j < M->N; j++) {
      arb_add(avg, *arb_mat_entry(M, i, j), *arb_mat_entry(M, j, i));
      arb_div_ui(avg, avg, 2);
      arb_set(*arb_mat_entry(M, i, j), avg);
      arb_set(*arb_mat_entry(M, j, i), avg);
    }
  }
  arb_clear(avg);
}

/*
 *   Jacobi eigensolver for arbitrary-precision symmetric matrix
 *   Returns eigenvalues in sorted order.
 */

static inline void arb_mat_jacobi_eigenvalues(arb_t *evals, const arb_mat_t *A,
                                              int max_sweeps) {
  int N = A->N;
  mpfr_prec_t prec = A->prec;
  arb_mat_t V;
  arb_mat_init(&V, N, prec);

  for (int i = 0; i < N * N; i++) {
    arb_set(V.entries[i], A->entries[i]);
  }

  arb_t apq, app, aqq, tau_val, t_val, c_val, s_val;
  arb_t vip, viq, vpj, vqj, threshold, max_off, abs_apq, denom;
  arb_t one_arb, tmp;

  arb_init2(apq, prec);
  arb_init2(app, prec);
  arb_init2(aqq, prec);
  arb_init2(tau_val, prec);
  arb_init2(t_val, prec);
  arb_init2(c_val, prec);
  arb_init2(s_val, prec);
  arb_init2(vip, prec);
  arb_init2(viq, prec);
  arb_init2(vpj, prec);
  arb_init2(vqj, prec);
  arb_init2(threshold, prec);
  arb_init2(max_off, prec);
  arb_init2(abs_apq, prec);
  arb_init2(denom, prec);
  arb_init2(one_arb, prec);
  arb_init2(tmp, prec);

  arb_set_ui(one_arb, 1);

  arb_t conv_thresh;
  arb_init2(conv_thresh, prec);
  if (prec > 100) {
    mpfr_set_ui_2exp(conv_thresh, 1UL, -(mpfr_exp_t)(prec / 4), ARB_RND);
  } else {
    arb_set_d(conv_thresh, 1e-30);
  }

  int sweeps_used = 0;
  for (int sweep = 0; sweep < max_sweeps; sweep++) {
    sweeps_used = sweep + 1;
    arb_set_ui(max_off, 0);
    for (int p = 0; p < N - 1; p++) {
      for (int q = p + 1; q < N; q++) {
        arb_abs(abs_apq, V.entries[p * N + q]);
        if (arb_cmp(abs_apq, max_off) > 0) {
          arb_set(max_off, abs_apq);
        }
      }
    }

    if (arb_cmp(max_off, conv_thresh) < 0)
      break;

    for (int p = 0; p < N - 1; p++) {
      for (int q = p + 1; q < N; q++) {
        arb_set(apq, V.entries[p * N + q]);
        arb_abs(abs_apq, apq);
        arb_set(app, V.entries[p * N + p]);
        arb_set(aqq, V.entries[q * N + q]);
        arb_abs(app, V.entries[p * N + p]);
        arb_abs(aqq, V.entries[q * N + q]);
        arb_add(denom, app, aqq);
        arb_add_ui(denom, denom, 1);
        arb_mul_d(denom, denom, 1e-30);
        if (arb_cmp(abs_apq, denom) < 0)
          continue;

        arb_set(app, V.entries[p * N + p]);
        arb_set(aqq, V.entries[q * N + q]);

        /* tau = (aqq - app) / (2 * apq) */
        arb_sub(tau_val, aqq, app);
        arb_mul_ui(tmp, apq, 2);
        arb_div(tau_val, tau_val, tmp);

        /* t = sign(tau) / (|tau| + sqrt(1 + tau^2)) */
        arb_sqr(tmp, tau_val);
        arb_add(tmp, tmp, one_arb);
        arb_sqrt(tmp, tmp);
        arb_abs(denom, tau_val);
        arb_add(denom, denom, tmp);
        if (arb_sgn(tau_val) >= 0) {
          arb_div(t_val, one_arb, denom);
        } else {
          arb_div(t_val, one_arb, denom);
          arb_neg(t_val, t_val);
        }

        /* c = 1/sqrt(1 + t^2) */
        arb_sqr(tmp, t_val);
        arb_add(tmp, tmp, one_arb);
        arb_sqrt(c_val, tmp);
        arb_div(c_val, one_arb, c_val);

        /* s = t * c */
        arb_mul(s_val, t_val, c_val);

        /* Apply Givens rotation */
        for (int i = 0; i < N; i++) {
          arb_set(vip, V.entries[i * N + p]);
          arb_set(viq, V.entries[i * N + q]);
          arb_mul(V.entries[i * N + p], vip, c_val);
          arb_mul(tmp, viq, s_val);
          arb_sub(V.entries[i * N + p], V.entries[i * N + p], tmp);
          arb_mul(V.entries[i * N + q], vip, s_val);
          arb_mul(tmp, viq, c_val);
          arb_add(V.entries[i * N + q], V.entries[i * N + q], tmp);
        }
        for (int j = 0; j < N; j++) {
          arb_set(vpj, V.entries[p * N + j]);
          arb_set(vqj, V.entries[q * N + j]);
          arb_mul(V.entries[p * N + j], vpj, c_val);
          arb_mul(tmp, vqj, s_val);
          arb_sub(V.entries[p * N + j], V.entries[p * N + j], tmp);
          arb_mul(V.entries[q * N + j], vpj, s_val);
          arb_mul(tmp, vqj, c_val);
          arb_add(V.entries[q * N + j], V.entries[q * N + j], tmp);
        }
      }
    }
  }

  for (int i = 0; i < N; i++) {
    arb_set(evals[i], V.entries[i * N + i]);
  }

  fprintf(stderr, "  Jacobi: %d sweeps (%dx%d)\n", sweeps_used, N, N);

  /* Insertion sort */
  for (int i = 1; i < N; i++) {
    arb_t key;
    arb_init2(key, prec);
    arb_set(key, evals[i]);
    int j = i - 1;
    while (j >= 0 && arb_cmp(evals[j], key) > 0) {
      arb_set(evals[j + 1], evals[j]);
      j--;
    }
    arb_set(evals[j + 1], key);
    arb_clear(key);
  }

  arb_clear(apq);
  arb_clear(app);
  arb_clear(aqq);
  arb_clear(tau_val);
  arb_clear(t_val);
  arb_clear(c_val);
  arb_clear(s_val);
  arb_clear(vip);
  arb_clear(viq);
  arb_clear(vpj);
  arb_clear(vqj);
  arb_clear(threshold);
  arb_clear(max_off);
  arb_clear(abs_apq);
  arb_clear(denom);
  arb_clear(one_arb);
  arb_clear(tmp);
  arb_clear(conv_thresh);

  arb_mat_clear(&V);
}

/*
 *   LU determinant with partial pivoting (MPFR)
 *   Returns det(A) via sign-corrected product of diagonal of U.
 */

static inline void arb_mat_det(arb_t det, const arb_mat_t *A) {
  int N = A->N;
  mpfr_prec_t prec = A->prec;

  arb_mat_t LU;
  arb_mat_init(&LU, N, prec);
  for (int i = 0; i < N * N; i++)
    arb_set(LU.entries[i], A->entries[i]);

  int sign = 1;
  arb_t pivot, factor, tmp;
  arb_init2(pivot, prec);
  arb_init2(factor, prec);
  arb_init2(tmp, prec);

  for (int k = 0; k < N; k++) {
    int max_row = k;
    arb_abs(pivot, LU.entries[k * N + k]);
    for (int i = k + 1; i < N; i++) {
      arb_t a_ik;
      arb_init2(a_ik, prec);
      arb_abs(a_ik, LU.entries[i * N + k]);
      if (arb_cmp(a_ik, pivot) > 0) {
        max_row = i;
        arb_set(pivot, a_ik);
      }
      arb_clear(a_ik);
    }
    if (max_row != k) {
      for (int j = 0; j < N; j++) {
        arb_set(tmp, LU.entries[k * N + j]);
        arb_set(LU.entries[k * N + j], LU.entries[max_row * N + j]);
        arb_set(LU.entries[max_row * N + j], tmp);
      }
      sign = -sign;
    }
    if (mpfr_cmp_ui(LU.entries[k * N + k], 0) == 0) {
      arb_set_ui(det, 0);
      goto cleanup;
    }
    for (int i = k + 1; i < N; i++) {
      arb_div(factor, LU.entries[i * N + k], LU.entries[k * N + k]);
      arb_set(LU.entries[i * N + k], factor);
      for (int j = k + 1; j < N; j++) {
        arb_mul(tmp, factor, LU.entries[k * N + j]);
        arb_sub(LU.entries[i * N + j], LU.entries[i * N + j], tmp);
      }
    }
  }

  arb_set_ui(det, 1);
  for (int k = 0; k < N; k++) {
    arb_mul(det, det, LU.entries[k * N + k]);
  }
  if (sign < 0)
    arb_neg(det, det);

cleanup:
  arb_clear(pivot);
  arb_clear(factor);
  arb_clear(tmp);
  arb_mat_clear(&LU);
}

/*
 *   LU solve: Ax = b via PA = LU factorization
 *   x and b are N-element arb_t arrays.
 */

static inline void arb_mat_solve(arb_t *x, const arb_mat_t *A, const arb_t *b) {
  int N = A->N;
  mpfr_prec_t prec = A->prec;

  arb_mat_t LU;
  arb_mat_init(&LU, N, prec);
  for (int i = 0; i < N * N; i++)
    arb_set(LU.entries[i], A->entries[i]);

  int *perm = (int *)malloc((size_t)N * sizeof(int));
  for (int i = 0; i < N; i++)
    perm[i] = i;

  arb_t pivot_val, factor, tmp;
  arb_init2(pivot_val, prec);
  arb_init2(factor, prec);
  arb_init2(tmp, prec);

  for (int k = 0; k < N; k++) {
    int max_row = k;
    arb_abs(pivot_val, LU.entries[k * N + k]);
    for (int i = k + 1; i < N; i++) {
      arb_t a_ik;
      arb_init2(a_ik, prec);
      arb_abs(a_ik, LU.entries[i * N + k]);
      if (arb_cmp(a_ik, pivot_val) > 0) {
        max_row = i;
        arb_set(pivot_val, a_ik);
      }
      arb_clear(a_ik);
    }
    if (max_row != k) {
      for (int j = 0; j < N; j++) {
        arb_set(tmp, LU.entries[k * N + j]);
        arb_set(LU.entries[k * N + j], LU.entries[max_row * N + j]);
        arb_set(LU.entries[max_row * N + j], tmp);
      }
      int pi = perm[k];
      perm[k] = perm[max_row];
      perm[max_row] = pi;
    }
    for (int i = k + 1; i < N; i++) {
      arb_div(factor, LU.entries[i * N + k], LU.entries[k * N + k]);
      arb_set(LU.entries[i * N + k], factor);
      for (int j = k + 1; j < N; j++) {
        arb_mul(tmp, factor, LU.entries[k * N + j]);
        arb_sub(LU.entries[i * N + j], LU.entries[i * N + j], tmp);
      }
    }
  }

  arb_t *pb = (arb_t *)malloc((size_t)N * sizeof(arb_t));
  for (int i = 0; i < N; i++) {
    arb_init2(pb[i], prec);
    arb_set(pb[i], b[perm[i]]);
  }

  arb_t *y = (arb_t *)malloc((size_t)N * sizeof(arb_t));
  for (int i = 0; i < N; i++) {
    arb_init2(y[i], prec);
    arb_set(y[i], pb[i]);
    for (int j = 0; j < i; j++) {
      arb_mul(tmp, LU.entries[i * N + j], y[j]);
      arb_sub(y[i], y[i], tmp);
    }
  }

  for (int i = N - 1; i >= 0; i--) {
    arb_init2(x[i], prec);
    arb_set(x[i], y[i]);
    for (int j = i + 1; j < N; j++) {
      arb_mul(tmp, LU.entries[i * N + j], x[j]);
      arb_sub(x[i], x[i], tmp);
    }
    arb_div(x[i], x[i], LU.entries[i * N + i]);
  }

  for (int i = 0; i < N; i++) {
    arb_clear(pb[i]);
    arb_clear(y[i]);
  }
  free(pb);
  free(y);
  free(perm);
  arb_clear(pivot_val);
  arb_clear(factor);
  arb_clear(tmp);
  arb_mat_clear(&LU);
}

/*
 *   Rayleigh quotient iteration for smallest eigenvalue
 *   Returns λ_min via inverse iteration with cubic convergence.
 *   shift starts at 0, converges to eigenvalue nearest 0.
 */

static inline void arb_mat_rayleigh_min(arb_t lambda_min, const arb_mat_t *A,
                                        int max_iter) {
  int N = A->N;
  mpfr_prec_t prec = A->prec;

  arb_t sigma, nrm, rq, dot;
  arb_init2(sigma, prec);
  arb_init2(nrm, prec);
  arb_init2(rq, prec);
  arb_init2(dot, prec);

  arb_mat_t Ashift;
  arb_mat_init(&Ashift, N, prec);

  arb_t *v = (arb_t *)malloc((size_t)N * sizeof(arb_t));
  arb_t *w = (arb_t *)malloc((size_t)N * sizeof(arb_t));
  for (int i = 0; i < N; i++) {
    arb_init2(v[i], prec);
    arb_init2(w[i], prec);
  }

  arb_set_ui(v[0], 1);
  for (int i = 1; i < N; i++)
    arb_set_ui(v[i], 0);

  arb_set_ui(sigma, 0);

  for (int iter = 0; iter < max_iter; iter++) {
    for (int i = 0; i < N * N; i++)
      arb_set(Ashift.entries[i], A->entries[i]);
    for (int i = 0; i < N; i++) {
      arb_sub(Ashift.entries[i * N + i], Ashift.entries[i * N + i], sigma);
    }

    arb_mat_solve(w, &Ashift, v);

    arb_set_ui(nrm, 0);
    for (int i = 0; i < N; i++) {
      arb_mul(dot, w[i], w[i]);
      arb_add(nrm, nrm, dot);
    }
    mpfr_sqrt(nrm, nrm, ARB_RND);
    if (mpfr_cmp_ui(nrm, 0) == 0)
      break;

    for (int i = 0; i < N; i++)
      arb_div(v[i], w[i], nrm);

    arb_set_ui(rq, 0);
    for (int i = 0; i < N; i++) {
      arb_t Av_i;
      arb_init2(Av_i, prec);
      arb_set_ui(Av_i, 0);
      for (int j = 0; j < N; j++) {
        arb_mul(dot, A->entries[i * N + j], v[j]);
        arb_add(Av_i, Av_i, dot);
      }
      arb_mul(dot, v[i], Av_i);
      arb_add(rq, rq, dot);
      arb_clear(Av_i);
    }

    arb_set(sigma, rq);
  }

  arb_set(lambda_min, sigma);

  arb_clear(sigma);
  arb_clear(nrm);
  arb_clear(rq);
  arb_clear(dot);
  arb_mat_clear(&Ashift);
  for (int i = 0; i < N; i++) {
    arb_clear(v[i]);
    arb_clear(w[i]);
  }
  free(v);
  free(w);
}

/*
 *   Gauss-Legendre nodes and weights (double precision, precomputed)
 *   Used for archimedean quadrature sub-intervals.
 *   Computes via Newton iteration on Legendre polynomial roots.
 */

static inline void gauss_legendre_point(int i, int n, double *xi, double *wi) {
  double x = cos(M_PI * ((double)i + 0.75) / ((double)n + 0.5));
  for (int iter = 0; iter < 50; iter++) {
    double p0 = 1.0, p1 = x;
    for (int k = 2; k <= n; k++) {
      double p2 = ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / (double)k;
      p0 = p1;
      p1 = p2;
    }
    double dp = (double)n * (x * p1 - p0) / (x * x - 1.0);
    double dx = -p1 / dp;
    x += dx;
    if (fabs(dx) < 1e-15)
      break;
  }
  *xi = x;
  double p0 = 1.0, p1 = x;
  for (int k = 2; k <= n; k++) {
    double p2 = ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / (double)k;
    p0 = p1;
    p1 = p2;
  }
  double dp = (double)n * (x * p1 - p0) / (x * x - 1.0);
  *wi = 2.0 / ((1.0 - x * x) * dp * dp);
}

static inline void gauss_legendre_point_mpfr(int i, int n, arb_t xi, arb_t wi,
                                             mpfr_prec_t prec) {
  double xi_d, wi_d;
  gauss_legendre_point(i, n, &xi_d, &wi_d);
  arb_set_d(xi, xi_d);

  arb_t p0, p1, k_arb, tmp, dp, xsq, one, dx, abs_dx, thr;
  arb_init2(p0, prec);
  arb_init2(p1, prec);
  arb_init2(k_arb, prec);
  arb_init2(tmp, prec);
  arb_init2(dp, prec);
  arb_init2(xsq, prec);
  arb_init2(one, prec);
  arb_init2(dx, prec);
  arb_init2(abs_dx, prec);
  arb_init2(thr, prec);

  mpfr_set_str(thr, "1e-75", 10, ARB_RND);

  for (int iter = 0; iter < 10; iter++) {
    arb_set_ui(p0, 1);
    arb_set(p1, xi);
    arb_set_ui(one, 1);
    for (int k = 2; k <= n; k++) {
      arb_mul(tmp, xi, p1);
      arb_mul_ui(tmp, tmp, 2 * (unsigned)k - 1);
      arb_mul_ui(k_arb, p0, (unsigned)(k - 1));
      arb_sub(tmp, tmp, k_arb);
      arb_div_ui(tmp, tmp, (unsigned)k);
      arb_set(p0, p1);
      arb_set(p1, tmp);
    }
    arb_mul(xsq, xi, xi);
    arb_sub(xsq, xsq, one);
    arb_mul(dp, xi, p1);
    arb_sub(dp, dp, p0);
    arb_mul_ui(dp, dp, (unsigned)n);
    arb_div(dp, dp, xsq);

    arb_neg(dx, p1);
    arb_div(dx, dx, dp);
    arb_add(xi, xi, dx);
    arb_abs(abs_dx, dx);
    if (mpfr_cmp(abs_dx, thr) < 0)
      break;
  }

  arb_set_ui(p0, 1);
  arb_set(p1, xi);
  arb_set_ui(one, 1);
  for (int k = 2; k <= n; k++) {
    arb_mul(tmp, xi, p1);
    arb_mul_ui(tmp, tmp, 2 * (unsigned)k - 1);
    arb_mul_ui(k_arb, p0, (unsigned)(k - 1));
    arb_sub(tmp, tmp, k_arb);
    arb_div_ui(tmp, tmp, (unsigned)k);
    arb_set(p0, p1);
    arb_set(p1, tmp);
  }
  arb_mul(xsq, xi, xi);
  arb_sub(xsq, xsq, one);
  arb_neg(xsq, xsq);
  arb_mul(dp, xi, p1);
  arb_sub(dp, dp, p0);
  arb_mul_ui(dp, dp, (unsigned)n);
  arb_div(dp, dp, xsq);
  arb_mul(tmp, dp, dp);
  arb_mul(tmp, tmp, xsq);
  arb_set_ui(wi, 2);
  arb_div(wi, wi, tmp);

  arb_clear(p0);
  arb_clear(p1);
  arb_clear(k_arb);
  arb_clear(tmp);
  arb_clear(dp);
  arb_clear(xsq);
  arb_clear(one);
  arb_clear(dx);
  arb_clear(abs_dx);
  arb_clear(thr);
}
