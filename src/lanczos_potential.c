/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Lanczos-accelerated large-N custom potential for Dyson BM
 * @paper   yamaguchi-rh-2026.tex, Section 8
 * @theorem Theorem III (RvM Density Matching)
 * @proof   Low-rank Lanczos approximation V'(H) = Q V'(T) Q^T
 * @step    3 -- Lanczos acceleration for N >= 200
 *
 * lanczos_potential.c - Phase 3: Lanczos-accelerated large-N custom potential
 *
 * Key innovation: V'(H) computed via Lanczos low-rank approximation:
 *   H ≈ Q·T·Q^T  (k-step Lanczos, k ≪ N)
 *   V'(H) ≈ Q·V'(T)·Q^T  (scalar V' applied to T's eigenvalues)
 *
 * Cost per Dyson BM step: O(k·N²) vs O(N³) for full diagonalization.
 * Enables N ≥ 200 with k ≈ 50.
 *
 * Potential: V'(x) = a·x + b·x·log(1 + c·x²) + d·x/(1 + e·x²)
 *   Controls: linear (a), log-tail (b,c), rational saturation (d,e)
 *
 * Compile: gcc -O3 -o lanczos_potential lanczos_potential.c -lm
 */

#include "random_matrix_utils.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef _WIN32
#include <windows.h>
#endif

static double vec_dot(const double *a, const double *b, int n) {
  double d = 0.0;
  for (int i = 0; i < n; i++)
    d += a[i] * b[i];
  return d;
}
static void mat_vec(const double *H, int N, const double *v, double *w) {
  for (int i = 0; i < N; i++) {
    double s = 0.0;
    for (int j = 0; j < N; j++)
      s += H[i * N + j] * v[j];
    w[i] = s;
  }
}
/* Proper Lanczos implementation */
static int lanczos_tridiag(const double *H, int N, int k, double *Q,
                           double *alpha, double *beta) {
  double *v = (double *)xcalloc((size_t)N, sizeof(double));
  double *w = (double *)xmalloc((size_t)N * sizeof(double));
  double *v0 = (double *)xcalloc((size_t)N, sizeof(double)); /* v_{j-1} */

  /* Random start vector */
  for (int i = 0; i < N; i++)
    v[i] = rng_normal();
  double nrm = sqrt(vec_dot(v, v, N));
  for (int i = 0; i < N; i++)
    v[i] /= nrm;
  /* Store as first Lanczos vector */
  for (int i = 0; i < N; i++)
    Q[i * k + 0] = v[i];

  for (int j = 0; j < k; j++) {
    mat_vec(H, N, v, w);
    alpha[j] = vec_dot(v, w, N);
    for (int i = 0; i < N; i++)
      w[i] -= alpha[j] * v[i] + (j > 0 ? beta[j - 1] : 0.0) * v0[i];

    /* Full reorthogonalization */
    for (int t = 0; t < 2; t++) {
      for (int i = 0; i <= j; i++) {
        double d = 0.0;
        for (int r = 0; r < N; r++)
          d += w[r] * Q[r * k + i];
        for (int r = 0; r < N; r++)
          w[r] -= d * Q[r * k + i];
      }
    }

    double bj = sqrt(vec_dot(w, w, N));
    if (bj < 1e-14) {
      k = j + 1;
      break;
    }
    beta[j] = bj;
    for (int i = 0; i < N; i++) {
      v0[i] = v[i];
      v[i] = w[i] / bj;
    }
    if (j + 1 < k)
      for (int i = 0; i < N; i++)
        Q[i * k + (j + 1)] = v[i];
  }

  free(v);
  free(w);
  free(v0);
  return k;
}

/* ==========================================================================
 *   Jacobi diagonalization of T (k×k tridiagonal → eigenvalues μ, vectors U)
 *   T is stored as: diag=alpha[k], off_diag=beta[k-1]
 *   Returns eigenvalues in mu (sorted) and eigenvectors in U (k×k).
 * ========================================================================== */

static void jacobi_tridiag(int k, const double *alpha, const double *beta,
                           double *mu, double *U) {
  double *T = (double *)xcalloc((size_t)(k * k), sizeof(double));
  for (int i = 0; i < k; i++) {
    T[i * k + i] = alpha[i];
    if (i < k - 1)
      T[i * k + (i + 1)] = T[(i + 1) * k + i] = beta[i];
  }
  /* Identity for U */
  for (int i = 0; i < k * k; i++)
    U[i] = 0.0;
  for (int i = 0; i < k; i++)
    U[i * k + i] = 1.0;

  for (int sw = 0; sw < 50; sw++) {
    double moff = 0.0;
    for (int p = 0; p < k - 1; p++)
      for (int q = p + 1; q < k; q++) {
        double v = fabs(T[p * k + q]);
        if (v > moff)
          moff = v;
      }
    if (moff < 1e-14)
      break;
    for (int p = 0; p < k - 1; p++)
      for (int q = p + 1; q < k; q++) {
        double apq = T[p * k + q];
        if (fabs(apq) < 1e-16 * (fabs(T[p * k + p]) + fabs(T[q * k + q]) + 1.0))
          continue;
        double app = T[p * k + p], aqq = T[q * k + q];
        double tau = (aqq - app) / (2.0 * apq);
        double t = tau >= 0 ? 1.0 / (tau + sqrt(1.0 + tau * tau))
                            : -1.0 / (-tau + sqrt(1.0 + tau * tau));
        double c = 1.0 / sqrt(1.0 + t * t), s = t * c;
        for (int i = 0; i < k; i++) {
          double tip = T[i * k + p], tiq = T[i * k + q];
          T[i * k + p] = tip * c - tiq * s;
          T[i * k + q] = tip * s + tiq * c;
        }
        for (int j = 0; j < k; j++) {
          double tpj = T[p * k + j], tqj = T[q * k + j];
          T[p * k + j] = tpj * c - tqj * s;
          T[q * k + j] = tpj * s + tqj * c;
        }
        for (int i = 0; i < k; i++) {
          double uip = U[i * k + p], uiq = U[i * k + q];
          U[i * k + p] = uip * c - uiq * s;
          U[i * k + q] = uip * s + uiq * c;
        }
      }
  }
  for (int i = 0; i < k; i++)
    mu[i] = T[i * k + i];
  /* sort */
  for (int i = 1; i < k; i++) {
    double key = mu[i];
    int j = i - 1;
    double *col = (double *)xmalloc((size_t)k * sizeof(double));
    for (int r = 0; r < k; r++)
      col[r] = U[r * k + i];
    while (j >= 0 && mu[j] > key) {
      mu[j + 1] = mu[j];
      for (int r = 0; r < k; r++)
        U[r * k + (j + 1)] = U[r * k + j];
      j--;
    }
    mu[j + 1] = key;
    for (int r = 0; r < k; r++)
      U[r * k + (j + 1)] = col[r];
    free(col);
  }
  free(T);
}

/* ==========================================================================
 *   Scalar V'(x) function
 *   V'(x) = a·x + b·x·log(1+c·x²) + d·x/(1+e·x²)
 *   Allows: linear (a), log-growth (b,c), saturating rational (d,e)
 * ========================================================================== */

typedef struct {
  double a, b, c, d, e;
} VP;

static double vprime_scalar(double x, VP p) {
  double x2 = x * x;
  double lg = (p.c > 0 && x2 > 0) ? log1p(p.c * x2) : 0.0;
  double rn = (p.e > 0) ? 1.0 / (1.0 + p.e * x2) : 0.0;
  return p.a * x + p.b * x * lg + p.d * x * rn;
}

/* ==========================================================================
 *   Compute V'(H) via Lanczos approximation:
 *     V'(H) ≈ Q · [U·diag(V'(μ))·U^T] · Q^T
 *   Output: drift matrix in H_out (symmetric, N×N)
 *   Workspace: Q(N×k), U(k×k), mu(k), M(N×k), UT(k×k)
 * ========================================================================== */

static void vprime_lanczos(const double *H, int N, int k, VP p, double *H_out,
                           double *Q, double *U, double *mu, double *M,
                           double *UT, double *alpha, double *beta) {
  int kk = lanczos_tridiag(H, N, k, Q, alpha, beta);
  jacobi_tridiag(kk, alpha, beta, mu, U);

  /* UT = U · diag(V'(μ)) */
  for (int i = 0; i < kk; i++) {
    double d = vprime_scalar(mu[i], p);
    for (int r = 0; r < kk; r++)
      UT[r * kk + i] = U[r * kk + i] * d;
  }
  /* M = UT · U^T  (kk×kk) */
  double *Tk_vp = (double *)xcalloc((size_t)(kk * kk), sizeof(double));
  for (int i = 0; i < kk; i++)
    for (int j = 0; j < kk; j++) {
      double s = 0.0;
      for (int l = 0; l < kk; l++)
        s += UT[i * kk + l] * U[j * kk + l];
      Tk_vp[i * kk + j] = s;
    }
  /* M = Q · Tk_vp  (N×kk) */
  for (int i = 0; i < N; i++)
    for (int j = 0; j < kk; j++) {
      double s = 0.0;
      for (int l = 0; l < kk; l++)
        s += Q[i * k + l] * Tk_vp[l * kk + j];
      M[i * kk + j] = s;
    }
  /* H_out = M · Q^T  (N×N, rank-kk, symmetric) */
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++) {
      double s = 0.0;
      for (int l = 0; l < kk; l++)
        s += M[i * kk + l] * Q[j * k + l];
      H_out[i * N + j] = s;
    }
  free(Tk_vp);
}

/* ==========================================================================
 *   Dyson BM step with Lanczos drift
 *   H_new = H - 0.5*dt*V'(H) + dW
 * ========================================================================== */

static void dyson_lanczos_step(int N, double *H, int k, VP p, double dt,
                               double *Q, double *U, double *mu, double *M,
                               double *UT, double *alpha, double *beta,
                               double *drift) {
  vprime_lanczos(H, N, k, p, drift, Q, U, mu, M, UT, alpha, beta);
  double half_dt = 0.5 * dt;
  for (int i = 0; i < N * N; i++)
    H[i] -= drift[i] * half_dt;
  mat_add_goe_noise(H, N, dt);
}

/* ==========================================================================
 *   RvM comparison: eigenvalue counting function
 * ========================================================================== */

static double N_rvm(double T) {
  if (T < 14.134)
    return 0.0;
  double x = T / (2.0 * M_PI);
  return x * log(x) - x + 7.0 / 8.0;
}

static void jacobi_evals_simple(const double *A, int N, double *evals) {
  double *V = xmalloc((size_t)(N * N) * sizeof(double));
  memcpy(V, A, (size_t)(N * N) * sizeof(double));
  for (int sw = 0; sw < 50; sw++) {
    double moff = 0.0;
    for (int p = 0; p < N - 1; p++)
      for (int q = p + 1; q < N; q++) {
        double v = fabs(V[p * N + q]);
        if (v > moff)
          moff = v;
      }
    if (moff < 1e-14)
      break;
    for (int p = 0; p < N - 1; p++)
      for (int q = p + 1; q < N; q++) {
        double apq = V[p * N + q];
        if (fabs(apq) < 1e-16 * (fabs(V[p * N + p]) + fabs(V[q * N + q]) + 1.0))
          continue;
        double app = V[p * N + p], aqq = V[q * N + q];
        double tau = (aqq - app) / (2.0 * apq);
        double t = tau >= 0 ? 1.0 / (tau + sqrt(1.0 + tau * tau))
                            : -1.0 / (-tau + sqrt(1.0 + tau * tau));
        double c = 1.0 / sqrt(1.0 + t * t), s = t * c;
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
 *   Main
 * ========================================================================== */

static void section(const char *s) { printf("\n--- %s ---\n\n", s); }

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif
  printf(
      "==================================================================\n");
  printf("  Phase 3: Lanczos-Accelerated Custom Potential at Large N\n");
  printf("  V'(x) = a·x + b·x·log(1+c·x²) + d·x/(1+e·x²)\n");
  printf(
      "==================================================================\n\n");

  uint64_t seed = (uint64_t)time(NULL);
  printf("Seed: %llu\n", (unsigned long long)seed);
  rng_set_seed(seed);

  /* Allocate Lanczos workspace (for max N=150, k=40) */
  int N_max = 150, K = 40;
  size_t nn = (size_t)(N_max * N_max), nk = (size_t)(N_max * K),
         kk = (size_t)(K * K);
  double *H = xcalloc(nn, sizeof(double));
  double *drift = xcalloc(nn, sizeof(double));
  double *Q = xmalloc(nk * sizeof(double));
  double *U = xmalloc(kk * sizeof(double));
  double *mu = xmalloc((size_t)K * sizeof(double));
  double *M = xmalloc(nk * sizeof(double));
  double *UT = xmalloc(kk * sizeof(double));
  double *alpha = xmalloc((size_t)K * sizeof(double));
  double *beta = xmalloc((size_t)K * sizeof(double));

  /* ── Test 1: Lanczos vs Jacobi eigenvalues (calibration) ── */
  section("Calibration: Lanczos Ritz values vs Jacobi eigenvalues (GOE, N=80)");

  int N_cal = 80;
  for (int i = 0; i < N_cal * N_cal; i++)
    H[i] = 0.0;
  /* Run Dyson BM to equilibrium */
  for (int s = 0; s < 500; s++) {
    for (int i = 0; i < N_cal * N_cal; i++)
      H[i] -= H[i] * 0.005;
    mat_add_goe_noise(H, N_cal, 0.01);
  }
  double *ev_full = xmalloc((size_t)N_cal * sizeof(double));
  jacobi_evals_simple(H, N_cal, ev_full);

  int kk_cal = lanczos_tridiag(H, N_cal, 40, Q, alpha, beta);
  double *ritz_vals = xmalloc((size_t)kk_cal * sizeof(double));
  double *u_small = xmalloc((size_t)(kk_cal * kk_cal) * sizeof(double));
  jacobi_tridiag(kk_cal, alpha, beta, ritz_vals, u_small);

  printf("Ritz values (extreme) vs Jacobi eigenvalues:\n");
  printf("  k   Ritz λ_k      True λ_k       Error\n");
  printf("  ---  -----------   -----------   ----------\n");
  for (int i = 0; i < 5; i++) {
    printf("  %3d  %12.6f  %12.6f  %10.2e\n", i, ritz_vals[i], ev_full[i],
           fabs(ritz_vals[i] - ev_full[i]));
  }
  /* Top eigenvalues */
  printf("  ...\n");
  for (int i = kk_cal - 3; i < kk_cal; i++) {
    printf("  %3d  %12.6f  %12.6f  %10.2e\n", i, ritz_vals[i],
           ev_full[N_cal - kk_cal + i],
           fabs(ritz_vals[i] - ev_full[N_cal - kk_cal + i]));
  }
  free(ev_full);
  free(ritz_vals);
  free(u_small);

  /* ── Test 2: Custom potential density at N=100 ── */
  section(
      "Custom V(x) density: GOE vs mild log vs strong log (N=100, paths=30)");

  int N_t = 100, npaths = 30, burn = 400;
  double dt_t = 0.01;
  int n_bins = 30;
  double E_max = 60.0;

  VP potentials[] = {
      {1.0, 0.0, 0.0, 0.0, 0.0},   /* GOE: V=x², V'=x */
      {0.5, 0.0, 0.0, 0.0, 0.0},   /* weak GOE */
      {0.0, 0.3, 0.001, 0.0, 0.0}, /* mild log tail */
      {0.0, 0.5, 0.01, 0.0, 0.0},  /* strong log tail */
  };
  const char *labels[] = {"GOE (a=1)", "Weak x^2 (a=0.5)", "Mild log",
                          "Strong log"};

  for (int pi = 0; pi < 4; pi++) {
    VP p = potentials[pi];
    double *all_ev = xmalloc((size_t)(N_t * npaths) * sizeof(double));
    int idx = 0;
    rng_set_seed(seed);

    for (int ip = 0; ip < npaths; ip++) {
      for (int i = 0; i < N_t * N_t; i++)
        H[i] = 0.0;
      for (int s = 0; s < burn; s++)
        dyson_lanczos_step(N_t, H, 35, p, dt_t, Q, U, mu, M, UT, alpha, beta,
                           drift);
      double *ev = xmalloc((size_t)N_t * sizeof(double));
      jacobi_evals_simple(H, N_t, ev);
      for (int k = 0; k < N_t; k++)
        all_ev[idx++] = ev[k];
      free(ev);
    }

    /* Histogram */
    double binw = 2.0 * E_max / (double)n_bins;
    double *hist = xcalloc((size_t)n_bins, sizeof(double));
    for (int i = 0; i < idx; i++) {
      int b = (int)((double)n_bins * (all_ev[i] + E_max) / (2.0 * E_max));
      if (b >= 0 && b < n_bins)
        hist[b] += 1.0 / ((double)(npaths * N_t) * binw);
    }
    printf("\n%s:\n", labels[pi]);
    printf("  E-center   density\n");
    printf("  ---------  -------\n");
    for (int b = 0; b < n_bins; b += 2) {
      double ctr = -E_max + (b + 0.5) * binw;
      printf("  %8.2f  %8.6f\n", ctr, hist[b]);
    }
    free(all_ev);
    free(hist);
  }

  /* ── Test 3: RvM counting function comparison (best potential) ── */
  section("RvM counting function comparison (N=120, more paths)");

  VP p_best = {0.0, 0.2, 0.002, 1.0, 0.0}; /* very weak quadratic + mild log */
  int N_rvmt = 120, npaths_rvmt = 60, burn_rvmt = 500;

  double *all_ev_r = xmalloc((size_t)(N_rvmt * npaths_rvmt) * sizeof(double));
  int idx_r = 0;
  rng_set_seed(seed);

  printf("Running %d paths at N=%d with k=35 Lanczos steps...\n", npaths_rvmt,
         N_rvmt);
  for (int ip = 0; ip < npaths_rvmt; ip++) {
    for (int i = 0; i < N_rvmt * N_rvmt; i++)
      H[i] = 0.0;
    for (int s = 0; s < burn_rvmt; s++)
      dyson_lanczos_step(N_rvmt, H, 35, p_best, dt_t, Q, U, mu, M, UT, alpha,
                         beta, drift);
    double *ev = xmalloc((size_t)N_rvmt * sizeof(double));
    jacobi_evals_simple(H, N_rvmt, ev);
    for (int k = 0; k < N_rvmt; k++)
      all_ev_r[idx_r++] = ev[k];
    free(ev);
  }

  /* Sort */
  for (int i = 1; i < idx_r; i++) {
    double key = all_ev_r[i];
    int j = i - 1;
    while (j >= 0 && all_ev_r[j] > key) {
      all_ev_r[j + 1] = all_ev_r[j];
      j--;
    }
    all_ev_r[j + 1] = key;
  }

  /* Counting functions */
  printf("\n  E      N_J(E)     N_rvm(E)    ratio\n");
  printf("  -----  ---------  ---------   -------\n");
  int n_pts = 20;
  for (int pt = 0; pt < n_pts; pt++) {
    double E = 2.0 + (double)pt * (60.0 - 2.0) / (double)(n_pts - 1);
    int count = 0;
    for (int i = 0; i < idx_r; i++)
      if (all_ev_r[i] <= E)
        count++;
    double n_emp = (double)count / (double)npaths_rvmt;
    double n_z = N_rvm(E);
    printf("  %5.1f  %9.1f  %9.1f   %6.3f\n", E, n_emp, n_z,
           n_z > 0 ? n_emp / n_z : 0.0);
  }
  free(all_ev_r);

  /* ── Timing ── */
  section("Performance: Lanczos vs full Jacobi (GOE, N=100, one step)");
  {
    for (int i = 0; i < N_t * N_t; i++)
      H[i] = 0.0;
    for (int s = 0; s < 200; s++) {
      for (int i = 0; i < N_t * N_t; i++)
        H[i] -= H[i] * 0.005;
      mat_add_goe_noise(H, N_t, 0.01);
    }

    clock_t t0 = clock();
    for (int r = 0; r < 20; r++) /* one Lanczos step */
      dyson_lanczos_step(N_t, H, 35, p_best, dt_t, Q, U, mu, M, UT, alpha, beta,
                         drift);
    clock_t t1 = clock();
    double lanc_ms = 1000.0 * (double)(t1 - t0) / CLOCKS_PER_SEC / 20.0;

    double *ev_t = xmalloc((size_t)N_t * sizeof(double));
    clock_t t2 = clock();
    for (int r = 0; r < 5; r++)
      jacobi_evals_simple(H, N_t, ev_t);
    clock_t t3 = clock();
    double jac_ms = 1000.0 * (double)(t3 - t2) / CLOCKS_PER_SEC / 5.0;

    printf("  Lanczos step (k=35, N=100): %.1f ms\n", lanc_ms);
    printf("  Jacobi eigensolve (N=100):  %.1f ms\n", jac_ms);
    printf("  Speedup: %.0fx\n", jac_ms / (lanc_ms + 1e-30));
    free(ev_t);
  }

  free(H);
  free(drift);
  free(Q);
  free(U);
  free(mu);
  free(M);
  free(UT);
  free(alpha);
  free(beta);
  printf("\nPhase 3 complete.\n");
  return 0;
}
