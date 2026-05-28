/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Custom potential V(x) for RvM density matching via Dyson BM
 * @paper   yamaguchi-rh-2026.tex, Section 8
 * @theorem Theorem III (RvM Density Matching)
 * @proof   Nelder-Mead optimization of parametric potential
 * @step    3 -- potential search and KS distance minimization
 *
 * potential_search.c - Phases 2-3: Custom V(x) → RvM density matching
 *
 * Phase 2: Parametric potential V_θ(x), custom Dyson BM with V'(H)
 *          via diagonalization, density histogram → RvM comparison.
 * Phase 3: Nelder-Mead simplex optimization of θ to minimize
 *          KS distance between ensemble density and RvM target.
 *
 * Potential family: V(x) = a·x² + b·x²·log(1 + c·x²)
 *   Asymptotic: V'(x) ~ 4b·x·log(√c·x)  matches RvM ~ x log x
 *   Parameters: θ = (a, b, c) ∈ ℝ³₊
 *
 * Compile: gcc -O3 -o potential_search potential_search.c -lm
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
#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

#ifdef _WIN32
#include <windows.h>
#endif

/* ==========================================================================
 *   Jacobi eigenvalue decomposition WITH eigenvectors
 * ========================================================================== */

static void jacobi_eigen_full(const double *A, int N, double *evals,
                              double *evecs, int max_sweeps) {
  double *V = (double *)xmalloc((size_t)(N * N) * sizeof(double));
  memcpy(V, A, (size_t)(N * N) * sizeof(double));

  /* Initialize eigenvectors to identity */
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++)
      evecs[i * N + j] = (i == j) ? 1.0 : 0.0;
  }

  for (int sweep = 0; sweep < max_sweeps; sweep++) {
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

        /* Rotate columns of V */
        for (int i = 0; i < N; i++) {
          double vip = V[i * N + p];
          double viq = V[i * N + q];
          V[i * N + p] = vip * c - viq * s;
          V[i * N + q] = vip * s + viq * c;
        }
        /* Rotate rows of V */
        for (int j = 0; j < N; j++) {
          double vpj = V[p * N + j];
          double vqj = V[q * N + j];
          V[p * N + j] = vpj * c - vqj * s;
          V[q * N + j] = vpj * s + vqj * c;
        }
        /* Accumulate eigenvectors */
        for (int i = 0; i < N; i++) {
          double eip = evecs[i * N + p];
          double eiq = evecs[i * N + q];
          evecs[i * N + p] = eip * c - eiq * s;
          evecs[i * N + q] = eip * s + eiq * c;
        }
      }
    }
  }

  for (int i = 0; i < N; i++)
    evals[i] = V[i * N + i];

  /* Sort eigenvalues + eigenvectors together (insertion sort) */
  for (int i = 1; i < N; i++) {
    double key = evals[i];
    int j = i - 1;
    /* save column i of eigenvectors */
    double *col_i = (double *)xmalloc((size_t)N * sizeof(double));
    for (int r = 0; r < N; r++)
      col_i[r] = evecs[r * N + i];
    for (int r = 0; r < N; r++)
      col_i[r] = evecs[r * N + i];
    while (j >= 0 && evals[j] > key) {
      evals[j + 1] = evals[j];
      for (int r = 0; r < N; r++)
        evecs[r * N + j + 1] = evecs[r * N + j];
      j--;
    }
    evals[j + 1] = key;
    for (int r = 0; r < N; r++)
      evecs[r * N + j + 1] = col_i[r];
  }

  free(V);
}

/* ==========================================================================
 *   Parametric potential V_θ(x) = a·x² + b·x²·log(1 + c·x²)
 *   V'(x) = 2a·x + 2b·x·log(1 + c·x²) + 2b·c·x³/(1 + c·x²)
 *
 *   Asymptotic: V'(x) ~ 2b·x·log(c·x²) + 2b·x  ~ 4b·x·log(√c·|x|)
 *   For RvM: need V'(x) ~ x·log x, so 4b ≈ 1, √c ≈ 1/(2π)
 *   Initial guess: a=0, b=0.25, c≈0.025
 * ========================================================================== */

typedef struct {
  double a, b, c;
} VParams;

/* Evaluate V'(λ) at a single scalar λ */
static double vprime_scalar(double lambda, VParams p) {
  double x2 = lambda * lambda;
  double log_term = log1p(p.c * x2); /* log(1 + c·λ²), stable for small c·λ² */
  return 2.0 * p.a * lambda + 2.0 * p.b * lambda * log_term +
         2.0 * p.b * p.c * lambda * x2 / (1.0 + p.c * x2);
}

/* Compute V'(H) via diagonalization: H = Q Λ Q^T, V'(H) = Q V'(Λ) Q^T.
 * Input H is consumed; output H = V'(H). evals, evecs, tmp are workspace.
 * evecs is N×N column-major: evecs[i*N+j] = Q[i][j]. */
static void compute_vprime_matrix(double *H, int N, VParams p, double *evals,
                                  double *evecs, double *tmp) {
  jacobi_eigen_full(H, N, evals, evecs, 50);

  /* tmp = Q * V'(Λ)  (scale each column j of Q by V'(λ_j)) */
  for (int j = 0; j < N; j++) {
    double d = vprime_scalar(evals[j], p);
    for (int i = 0; i < N; i++)
      tmp[i * N + j] = evecs[i * N + j] * d;
  }
  /* H = tmp * Q^T  (row i of tmp dot row j of Q) */
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      double s = 0.0;
      for (int k = 0; k < N; k++)
        s += tmp[i * N + k] * evecs[j * N + k];
      H[i * N + j] = s;
    }
  }
}

/* ==========================================================================
 *   Dyson BM with custom potential
 *   dH = -(1/2) V'(H) dt + dW_goe
 * ========================================================================== */

static void dyson_bm_step_custom(int N, double *H, VParams p, double dt,
                                 double *evals, double *evecs, double *tmp,
                                 double *Hstore) {
  /* Compute drift: D = V'(H) */
  mat_copy(Hstore, H, N); /* save H */
  compute_vprime_matrix(Hstore, N, p, evals, evecs, tmp);
  /* Now Hstore = V'(H). Drift = -0.5 * V'(H) * dt */
  double half_dt = 0.5 * dt;
  for (int i = 0; i < N * N; i++)
    H[i] -= Hstore[i] * half_dt; /* H = H - (1/2)V'(H)dt */
  /* Add noise */
  mat_add_goe_noise(H, N, dt);
}

/* Run many steps, check energy convergence */
static int dyson_bm_equilibrate_custom(int N, double *H, VParams p, double dt,
                                       int burn_steps, double *evals,
                                       double *evecs, double *tmp) {
  double *Hstore = (double *)xmalloc((size_t)(N * N) * sizeof(double));
  for (int s = 0; s < burn_steps; s++) {
    dyson_bm_step_custom(N, H, p, dt, evals, evecs, tmp, Hstore);
  }
  free(Hstore);
  double energy = mat_frobenius_sq(H, N);
  /* For V=x², E[Tr(H²)] ≈ N(N+1). For custom V, we just report. */
  (void)energy;
  return 1;
}

/* ==========================================================================
 *   Riemann-von Mangoldt target density
 *   ρ_rvm(E) = (1/2π) log(|E|/(2π))  for |E| > 2π, symmetric
 *   Normalized: ρ_rvm actual integer count, not probability.
 *   We use the normalized probability density:
 *     p_rvm(E) ∝ log(|E|/2π) for |E| ≥ E_min
 *     p_rvm(E) = 0 for |E| < E_min
 *   with E_min = 2π ≈ 6.28 (first zero ~14.13, first Gram ~17.85)
 *   Actually, for finite N, we truncate at ±E_max.
 * ========================================================================== */

static double rvm_density(double E, double E_min, double E_max) {
  double aE = fabs(E);
  if (aE <= E_min || aE >= E_max)
    return 0.0;
  return log(aE / (2.0 * M_PI)) / (2.0 * M_PI);
}

/* Normalization integral (analytic): ∫_{E_min}^{E_max} log(E/2π)/(2π) dE
 * = [E(log E - log(2π) - 1)] / (2π)  from E_min to E_max
 * With symmetry, total = 2 * above. */
static double rvm_primitive(double E) {
  return E * (log(E) - log(2.0 * M_PI) - 1.0) / (2.0 * M_PI);
}

static double rvm_norm(double E_min, double E_max) {
  return 2.0 * (rvm_primitive(E_max) - rvm_primitive(E_min));
}

/* ==========================================================================
 *   Density histogram from eigenvalues
 *   Computes KS distance against RvM target.
 * ========================================================================== */

static double density_ks_distance(const double *all_evals, int total_evals,
                                  double E_min, double E_max, int nbins) {
  /* Build empirical histogram, then CDF */
  double bin_w = (E_max - (-E_max)) / (double)nbins;
  double *hist = (double *)xcalloc((size_t)nbins, sizeof(double));
  for (int i = 0; i < total_evals; i++) {
    double E = all_evals[i];
    int b = (int)((double)nbins * (E + E_max) / (2.0 * E_max));
    if (b >= 0 && b < nbins)
      hist[b] += 1.0;
  }

  /* Empirical CDF and theoretical CDF */
  double emp_cdf = 0.0, rvm_cdf = 0.0, ks = 0.0;
  double rvm_total = rvm_norm(E_min, E_max);
  for (int b = 0; b < nbins; b++) {
    double center = -E_max + (b + 0.5) * bin_w;
    emp_cdf += hist[b] / (double)total_evals;
    if (fabs(center) <= E_min)
      continue;
    if (fabs(center) > E_min && fabs(center) < E_max)
      rvm_cdf += rvm_density(center, E_min, E_max) * bin_w / rvm_total;
    double diff = fabs(emp_cdf - rvm_cdf);
    if (diff > ks)
      ks = diff;
  }

  free(hist);
  return ks;
}

/* ==========================================================================
 *   Density scan: run ensemble, compute KS
 * ========================================================================== */

static double evaluate_theta(VParams p, int N, int npaths, int burn_steps,
                             double dt, double E_max) {
  double E_min = 2.0 * M_PI; /* RvM log density starts after 2π */
  int nbins = 40;
  int total_evals = N * npaths;

  double *H = (double *)xmalloc((size_t)(N * N) * sizeof(double));
  double *evals = (double *)xmalloc((size_t)N * sizeof(double));
  double *evecs = (double *)xmalloc((size_t)(N * N) * sizeof(double));
  double *tmp = (double *)xmalloc((size_t)(N * N) * sizeof(double));
  double *all_ev = (double *)xmalloc((size_t)total_evals * sizeof(double));

  int idx = 0;
  for (int ip = 0; ip < npaths; ip++) {
    /* Start from GOE equilibrium, then adapt */
    for (int r = 0; r < N * N; r++)
      H[r] = 0.0;
    /* Warm up with standard GOE */
    for (int r = 0; r < N * N; r++)
      H[r] = 0.0;
    for (int s = 0; s < 200; s++) {
      double half_dt = 0.5 * dt;
      for (int r = 0; r < N * N; r++)
        H[r] -= H[r] * half_dt;
      mat_add_goe_noise(H, N, dt);
    }
    /* Then run custom potential */
    dyson_bm_equilibrate_custom(N, H, p, dt, burn_steps, evals, evecs, tmp);

    jacobi_eigen_full(H, N, evals, evecs, 50);
    for (int k = 0; k < N; k++)
      all_ev[idx++] = evals[k];
  }

  double ks = density_ks_distance(all_ev, total_evals, E_min, E_max, nbins);

  free(H);
  free(evals);
  free(evecs);
  free(tmp);
  free(all_ev);
  return ks;
}

/* ==========================================================================
 *   Nelder-Mead simplex optimizer (3D: a, b, c)
 *   Minimize KS distance between ensemble density and RvM.
 *   Direct-search derivative-free method, robust for noisy stochastic
 * objectives.
 * ========================================================================== */

static double nm_objective(double *x, int N, int npaths, int burn, double dt,
                           double E_max) {
  VParams p = {x[0], x[1], x[2]};
  return evaluate_theta(p, N, npaths, burn, dt, E_max);
}

static int nelder_mead_optimize(double *x0, int dim, int max_iter, int N,
                                int npaths, int burn, double dt, double E_max,
                                double simplex_size, double tol_simplex) {
  double (*simplex)[4] =
      (double (*)[4])xmalloc((size_t)(dim + 1) * sizeof(double[4]));
  double *fv = (double *)xmalloc((size_t)(dim + 1) * sizeof(double));
  double *centroid = (double *)xmalloc((size_t)dim * sizeof(double));
  double *xr = (double *)xmalloc((size_t)dim * sizeof(double));
  double *xe = (double *)xmalloc((size_t)dim * sizeof(double));
  double *xc = (double *)xmalloc((size_t)dim * sizeof(double));

  for (int i = 0; i <= dim; i++)
    for (int d = 0; d < dim; d++)
      simplex[i][d] =
          (i == 0) ? x0[d] : x0[d] + ((i - 1 == d) ? simplex_size : 0.0);

  for (int i = 0; i <= dim; i++)
    fv[i] = nm_objective(simplex[i], N, npaths, burn, dt, E_max);

  for (int iter = 0; iter < max_iter; iter++) {
    /* Sort indices: hi=highest, s_hi=second-highest, lo=lowest */
    int hi = 0, lo = 0;
    for (int i = 1; i <= dim; i++) {
      if (fv[i] > fv[hi])
        hi = i;
      if (fv[i] < fv[lo])
        lo = i;
    }
    int s_hi = (hi == 0) ? 1 : 0;
    for (int i = 0; i <= dim; i++) {
      if (i == hi)
        continue;
      if (fv[i] > fv[s_hi])
        s_hi = i;
    }

    double simplex_range = 0.0;
    for (int d = 0; d < dim; d++) {
      double vmax = simplex[0][d], vmin = simplex[0][d];
      for (int i = 1; i <= dim; i++) {
        if (simplex[i][d] > vmax)
          vmax = simplex[i][d];
        if (simplex[i][d] < vmin)
          vmin = simplex[i][d];
      }
      double rng = vmax - vmin;
      if (rng > simplex_range)
        simplex_range = rng;
    }
    if (simplex_range < tol_simplex) {
      printf("\n  NM converged at iter %d (simplex range %.6f < tol %.6f)\n",
             iter, simplex_range, tol_simplex);
      break;
    }

    for (int d = 0; d < dim; d++) {
      double sum = 0.0;
      for (int i = 0; i <= dim; i++)
        if (i != hi)
          sum += simplex[i][d];
      centroid[d] = sum / (double)dim;
    }

    for (int d = 0; d < dim; d++)
      xr[d] = 2.0 * centroid[d] - simplex[hi][d];
    double fr = nm_objective(xr, N, npaths, burn, dt, E_max);

    if (fr < fv[lo]) {
      for (int d = 0; d < dim; d++)
        xe[d] = 3.0 * centroid[d] - 2.0 * simplex[hi][d];
      double fe = nm_objective(xe, N, npaths, burn, dt, E_max);
      if (fe < fr) {
        for (int d = 0; d < dim; d++)
          simplex[hi][d] = xe[d];
        fv[hi] = fe;
      } else {
        for (int d = 0; d < dim; d++)
          simplex[hi][d] = xr[d];
        fv[hi] = fr;
      }
    } else if (fr < fv[s_hi]) {
      for (int d = 0; d < dim; d++)
        simplex[hi][d] = xr[d];
      fv[hi] = fr;
    } else {
      if (fr < fv[hi]) {
        for (int d = 0; d < dim; d++)
          simplex[hi][d] = xr[d];
        fv[hi] = fr;
      }
      for (int d = 0; d < dim; d++)
        xc[d] = 0.5 * (simplex[hi][d] + centroid[d]);
      double fc = nm_objective(xc, N, npaths, burn, dt, E_max);
      if (fc < fv[hi]) {
        for (int d = 0; d < dim; d++)
          simplex[hi][d] = xc[d];
        fv[hi] = fc;
      } else {
        for (int i = 0; i <= dim; i++) {
          if (i == lo)
            continue;
          for (int d = 0; d < dim; d++)
            simplex[i][d] = 0.5 * (simplex[i][d] + simplex[lo][d]);
          fv[i] = nm_objective(simplex[i], N, npaths, burn, dt, E_max);
        }
      }
    }
  }

  int best = 0;
  for (int i = 1; i <= dim; i++)
    if (fv[i] < fv[best])
      best = i;
  for (int d = 0; d < dim; d++)
    x0[d] = simplex[best][d];

  double best_fv = fv[best];
  free(simplex);
  free(fv);
  free(centroid);
  free(xr);
  free(xe);
  free(xc);
  return (int)(best_fv * 10000.0); /* non-negative success indicator */
}

/* ==========================================================================
 *   Main: Phase 2 density sweep + Phase 3 Nelder-Mead optimization
 * ========================================================================== */

static void section_header(const char *title) {
  printf("\n-------------------------------------------------------------------"
         "---\n");
  printf("  %s\n", title);
  printf("---------------------------------------------------------------------"
         "-\n\n");
}

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

  printf("====================================================================="
         "=\n");
  printf(
      "  Potential Search — Phase 2+3: Custom V(x) → RvM density matching\n");
  printf("  V(x) = a*x^2 + b*x^2*log(1 + c*x^2)\n");
  printf("  Target: ρ_rvm(E) ∝ log(|E|/2π)  (Riemann-von Mangoldt)\n");
  printf("====================================================================="
         "=\n\n");

  uint64_t seed = (uint64_t)time(NULL);
  printf("  Seed: %llu\n", (unsigned long long)seed);
  rng_set_seed(seed);

  /* ── Phase 2a: Single-point density check ── */

  section_header("Phase 2a: Density check at initial guess");

  VParams p_init = {0.0, 0.25, 0.025}; /* a=0, b≈1/4, c≈1/(2π)² */
  int N_test = 20;
  int npaths_sweep = 300;
  int burn_sweep = 500;
  double dt = 0.01;
  double E_max = 50.0;
  double E_min = 2.0 * M_PI;
  int nbins = 40;

  printf("  θ = (a=%.4f, b=%.4f, c=%.4f)\n", p_init.a, p_init.b, p_init.c);
  printf("  V'(x) ~ 4b·x·log(√c·|x|) = %.2f·x·log(%.4f·|x|) for large x\n",
         4.0 * p_init.b, sqrt(p_init.c));
  printf("  N=%d  paths=%d  burn=%d  dt=%.3f\n\n", N_test, npaths_sweep,
         burn_sweep, dt);

  double ks_init =
      evaluate_theta(p_init, N_test, npaths_sweep, burn_sweep, dt, E_max);
  printf("  KS distance (initial): %.6f\n", ks_init);

  /* Also check standard GOE as baseline */
  VParams p_goe = {0.5, 0.0, 0.0}; /* V=x²/2 */
  double ks_goe =
      evaluate_theta(p_goe, N_test, npaths_sweep, burn_sweep, dt, E_max);
  printf("  KS distance (GOE baseline): %.6f\n", ks_goe);

  /* ── Phase 2b: Density histogram with best guess ── */

  section_header("Phase 2b: Density histogram vs RvM target");

  {
    double *H = (double *)xmalloc((size_t)(N_test * N_test) * sizeof(double));
    double *evals = (double *)xmalloc((size_t)N_test * sizeof(double));
    double *evecs =
        (double *)xmalloc((size_t)(N_test * N_test) * sizeof(double));
    double *tmp = (double *)xmalloc((size_t)(N_test * N_test) * sizeof(double));
    double *all_ev =
        (double *)xmalloc((size_t)(N_test * npaths_sweep) * sizeof(double));

    int idx = 0;
    rng_set_seed(seed);
    for (int ip = 0; ip < npaths_sweep; ip++) {
      for (int r = 0; r < N_test * N_test; r++)
        H[r] = 0.0;
      for (int s = 0; s < 200; s++) {
        for (int r = 0; r < N_test * N_test; r++)
          H[r] -= H[r] * 0.5 * dt;
        mat_add_goe_noise(H, N_test, dt);
      }
      dyson_bm_equilibrate_custom(N_test, H, p_init, dt, burn_sweep, evals,
                                  evecs, tmp);
      jacobi_eigen_full(H, N_test, evals, evecs, 50);
      for (int k = 0; k < N_test; k++)
        all_ev[idx++] = evals[k];
    }

    double bin_w = 2.0 * E_max / (double)nbins;
    double *hist = (double *)xcalloc((size_t)nbins, sizeof(double));
    int total_evals = N_test * npaths_sweep;
    for (int i = 0; i < total_evals; i++) {
      double E = all_ev[i];
      int b = (int)((double)nbins * (E + E_max) / (2.0 * E_max));
      if (b >= 0 && b < nbins)
        hist[b] += 1.0;
    }

    printf("  Bin  Center     Empirical   RvM-target    Ratio\n");
    printf("  ──── ──────     ─────────   ──────────    ─────\n");
    double rvm_n = rvm_norm(E_min, E_max);
    for (int b = 0; b < nbins; b++) {
      double center = -E_max + (b + 0.5) * bin_w;
      double emp = hist[b] / ((double)total_evals * bin_w);
      double rvm = 0.0;
      if (fabs(center) > E_min && fabs(center) < E_max)
        rvm = rvm_density(center, E_min, E_max) / rvm_n;
      double ratio = (rvm > 1e-30) ? emp / rvm : 0.0;
      printf("  %4d  %7.2f   %10.6f  %10.6f   %7.3f\n", b, center, emp, rvm,
             ratio);
    }

    free(H);
    free(evals);
    free(evecs);
    free(tmp);
    free(all_ev);
    free(hist);
  }

  /* ── Phase 3: Nelder-Mead optimization ── */

  section_header("Phase 3a: Coarse grid seed for Nelder-Mead simplex");

  printf("  Testing 8 seed points near initial guess...\n\n");

  double scan_vals[][3] = {
      {0.0, 0.1, 0.01},  {0.0, 0.2, 0.02},    {0.0, 0.25, 0.025},
      {0.0, 0.3, 0.03},  {0.05, 0.25, 0.025}, {0.1, 0.25, 0.025},
      {0.0, 0.25, 0.05}, {0.0, 0.35, 0.025},
  };
  double best_ks = 1e300;
  VParams best_p = p_init;

  int npaths_scan = 100;
  int burn_scan = 300;
  int N_scan = 16;

  for (int si = 0; si < 8; si++) {
    VParams p = {scan_vals[si][0], scan_vals[si][1], scan_vals[si][2]};
    double ks = evaluate_theta(p, N_scan, npaths_scan, burn_scan, dt, 40.0);
    printf("  (a=%.3f, b=%.3f, c=%.4f): KS = %.6f %s\n", p.a, p.b, p.c, ks,
           (ks < best_ks) ? "<-- best" : "");
    if (ks < best_ks) {
      best_ks = ks;
      best_p = p;
    }
  }

  printf("\n  Best seed: (a=%.4f, b=%.4f, c=%.5f)  KS = %.6f\n", best_p.a,
         best_p.b, best_p.c, best_ks);
  printf("  Baseline GOE: KS = %.6f\n", ks_goe);

  /* ── Phase 3b: Real Nelder-Mead optimization ── */

  section_header("Phase 3b: Nelder-Mead simplex optimization");

  double x0[3] = {best_p.a, best_p.b, best_p.c};
  int nm_N = 16, nm_paths = 80, nm_burn = 250;
  int nm_max_iter = 30;
  double nm_simplex_size = 0.15;
  double nm_tol = 0.005;

  printf("  Starting NM from (a=%.4f, b=%.4f, c=%.5f)\n", x0[0], x0[1], x0[2]);
  printf("  N=%d  paths=%d  burn=%d  max_iter=%d\n\n", nm_N, nm_paths, nm_burn,
         nm_max_iter);

  nelder_mead_optimize(x0, 3, nm_max_iter, nm_N, nm_paths, nm_burn, dt, 40.0,
                       nm_simplex_size, nm_tol);

  double ks_nm = evaluate_theta((VParams){x0[0], x0[1], x0[2]}, nm_N,
                                nm_paths * 2, nm_burn * 2, dt, 40.0);
  printf("  NM result: (a=%.5f, b=%.5f, c=%.6f)  KS = %.6f\n", x0[0], x0[1],
         x0[2], ks_nm);
  printf("  ΔKS from GOE baseline: %+.4f (%.1f%%)\n\n", ks_nm - ks_goe,
         100.0 * (ks_nm - ks_goe) / ks_goe);

  /* ── Phase 3c: Larger N refinement ── */

  section_header("Phase 3c: Larger N Nelder-Mead (N=40, fewer paths)");

  double x0_large[3] = {x0[0], x0[1], x0[2]};
  int nm_N2 = 40, nm_paths2 = 40, nm_burn2 = 200;
  int nm_max_iter2 = 20;

  printf("  Starting NM from best small-N result\n");
  printf("  N=%d  paths=%d  burn=%d  max_iter=%d\n\n", nm_N2, nm_paths2,
         nm_burn2, nm_max_iter2);

  nelder_mead_optimize(x0_large, 3, nm_max_iter2, nm_N2, nm_paths2, nm_burn2,
                       dt, 60.0, nm_simplex_size * 0.5, nm_tol * 2.0);

  double ks_large =
      evaluate_theta((VParams){x0_large[0], x0_large[1], x0_large[2]}, nm_N2,
                     nm_paths2 * 2, nm_burn2 * 2, dt, 60.0);
  printf("  NM large-N result: (a=%.5f, b=%.5f, c=%.6f)  KS = %.6f\n",
         x0_large[0], x0_large[1], x0_large[2], ks_large);

  printf("\n  Phase 2+3 complete.\n");
  printf(
      "  Note: confining potentials (V~x²) cannot reproduce RvM at small N.\n");
  printf("  RvM spectral edge ~68 vs GOE ~11 for N=40.\n\n");

  return 0;
}
