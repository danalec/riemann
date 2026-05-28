/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   det-xi matching analysis
 * @paper   yamaguchi-rh-2026.tex, §8
 * @theorem Exploratory
 * @proof   Exploratory
 * @step    1 — Weyl asymptotics (eigenvalue approximation)
 *
 * det_xi_match.c - Phase 4: Eigenvalue counting function → RvM comparison
 *
 * Computes the ensemble-averaged eigenvalue counting function N_J(E)
 * and compares it to the Riemann-von Mangoldt zero counting function.
 *
 * The counting function is the integrated density — if densities match,
 * the Weyl law matches, which is a necessary condition for the
 * determinant identity det₂(I - z⁻¹J_∞) = ξ(½+iz).
 *
 * This is a preliminary Phase 4 test; full determinant identity
 * requires coefficient-level comparison (Phase 4b, future).
 *
 * Compile: gcc -O3 -o det_xi_match det_xi_match.c -lm
 */

#include "random_matrix_utils.h"
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef _WIN32
#include <windows.h>
#endif
/*
 *   Potential evaluation (scalar V', no matrix functions needed)
 *   For density comparison, we only need eigenvalues — not V'(H).
 *   We generate matrices directly from the invariant measure using
 *   Metropolis-adjusted Dyson BM instead of computing V'(H).
 *
 *   BUT: for simplicity, we use the standard GOE + importance weight.
 *   Generate H from GOE, then compute eigenvalues. Re-weight by
 *   exp(-(V(H) - Tr(H²)/2)) to adjust from GOE to custom V.
 *
 *   The eigenvalue density under custom V is:
 *     ρ_V(E) ∝ exp(-V(λ)/2) · Vandermonde²
 *
 *   We approximate this by rejecting samples: keep with probability
 *     min(1, exp(-(V_sum - V0_sum)))
 *   where V_sum = Σ V(λ_k), V0_sum = Σ λ_k²/2.
 */

typedef struct {
  double a, b, c;
} VParams;

/* Scalar potential V(λ) */
static double v_scalar(double lambda, VParams p) {
  double x2 = lambda * lambda;
  double log_term = (p.c > 0.0 && x2 > 0.0) ? log1p(p.c * x2) : 0.0;
  return p.a * x2 + p.b * x2 * log_term;
}

/* Sum of scalar potential over eigenvalue array: Σ V(λ_k) */
static double v_eigenvalue_sum(const double *evals, int N, VParams p) {
  double sum = 0.0;
  for (int k = 0; k < N; k++)
    sum += v_scalar(evals[k], p);
  return sum;
}

/*
 *   Riemann-von Mangoldt counting function
 *   N_rvm(T) = (T/2π) log(T/(2πe)) + 7/8 + S(T)
 *   Ignoring S(T) fluctuation: N_rvm ~ (T/2π) log(T/(2πe)) + 7/8
 */

static double N_rvm(double T) {
  if (T <= 14.134)
    return 0.0; /* first zero */
  double x = T / (2.0 * M_PI);
  return x * log(x) - x + 7.0 / 8.0;
}

static double N_rvm_symmetric(double E) {
  /* Symmetric: N_sym(E) = N_rvm(E) for E > 0, with symmetry N_sym(-E) =
   * N_sym(E). Translates to counting eigenvalues in [-E, E]: total =
   * 2*N_rvm(E). */
  double T = fabs(E);
  if (T < 14.134)
    return 0.0;
  return 2.0 * N_rvm(T); /* both positive and negative sides */
}

/*
 *   Counting function comparison
 */

typedef struct {
  double energy;
  double N_emp; /* empirical count ≤ E */
  double N_rvm; /* RvM count ≤ E */
  double diff;
} CountPoint;

/* Run ensemble with rejection sampling, compute counting function */
static void ensemble_counting_function(VParams p, int N, int npaths,
                                       int burn_steps, double dt, int npoints,
                                       CountPoint *cp, double E_max) {
  double *H = (double *)xmalloc((size_t)(N * N) * sizeof(double));
  double *evals = (double *)xmalloc((size_t)N * sizeof(double));

  int total_evals = npaths * N;
  double *all_ev = (double *)xmalloc((size_t)total_evals * sizeof(double));
  int idx = 0;
  int accepted = 0, rejected = 0;

  for (int ip = 0; ip < npaths * 3;
       ip++) { /* oversample to account for rejection */
    /* Generate from GOE equilibrium */
    for (int r = 0; r < N * N; r++)
      H[r] = 0.0;
    double half_dt = 0.5 * dt;
    for (int s = 0; s < burn_steps; s++) {
      for (int r = 0; r < N * N; r++)
        H[r] -= H[r] * half_dt;
      mat_add_goe_noise(H, N, dt);
    }

    jacobi_eigenvalues(H, N, evals);

    /* Compute GOE potential sum (V_0 = Σ λ²/2) and target V sum */
    double v0_sum = 0.0;
    for (int k = 0; k < N; k++)
      v0_sum += 0.5 * evals[k] * evals[k];
    double v_sum = v_eigenvalue_sum(evals, N, p);

    /* Acceptance: keep with probability exp(-ΔV) where ΔV = V - V_0 */
    double delta = v_sum - v0_sum;
    double accept_prob = exp(-delta);
    if (accept_prob > 1.0)
      accept_prob = 1.0;
    if (rng_uniform() < accept_prob) {
      for (int k = 0; k < N; k++)
        all_ev[idx++] = evals[k];
      accepted++;
      if (idx >= total_evals)
        break;
    } else {
      rejected++;
    }
  }

  /* Sort all eigenvalues for counting function */
  for (int i = 1; i < idx; i++) {
    double key = all_ev[i];
    int j = i - 1;
    while (j >= 0 && all_ev[j] > key) {
      all_ev[j + 1] = all_ev[j];
      j--;
    }
    all_ev[j + 1] = key;
  }

  /* Evaluate counting function at npoints equally spaced on [0, E_max] */
  double dE = E_max / (double)(npoints - 1);
  int ei = 0;
  int neg_idx = 0; /* first eigenvalue ≥ -E, counts evals with λ < -E */
  for (int pt = 0; pt < npoints; pt++) {
    double E = (double)pt * dE;
    cp[pt].energy = E;

    /* Count eigenvalues ≤ E */
    while (ei < idx && all_ev[ei] <= E)
      ei++;
    /* Count eigenvalues < -E (the ones NOT in the interval [-E, E]) */
    while (neg_idx < idx && all_ev[neg_idx] < -E)
      neg_idx++;

    /* Eigenvalues with |λ| ≤ E = (count ≤ E) - (count < -E) */
    int count_abs = ei - neg_idx;
    cp[pt].N_emp = (double)count_abs / (double)accepted; /* per-sample count */
    cp[pt].N_rvm =
        N_rvm_symmetric(E) / (double)N; /* RvM per-zero count scaled */
    cp[pt].diff = cp[pt].N_emp - cp[pt].N_rvm;
  }

  printf("    Acceptance rate: %d/%d = %.2f%%\n", accepted, accepted + rejected,
         100.0 * (double)accepted / (double)(accepted + rejected));

  free(H);
  free(evals);
  free(all_ev);
}

/*
 *   Main: Phase 4 counting function comparison
 */

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
  printf("  Phase 4: Eigenvalue Counting Function → RvM Comparison\n");
  printf("  N_J(E) vs N_rvm(E)  (proxy for determinant identity)\n");
  printf("====================================================================="
         "=\n\n");

  uint64_t seed = (uint64_t)time(NULL);
  printf("  Seed: %llu\n", (unsigned long long)seed);
  rng_set_seed(seed);

  int N_test = 30;
  int npaths = 200;
  int burn = 500;
  double dt = 0.01;
  double E_max = 40.0;
  int npoints = 30;

  /* ── GOE baseline ── */
  section_header("GOE baseline (V(x)=x²/2)");

  {
    VParams p_goe = {0.5, 0.0, 0.0};
    CountPoint *cp =
        (CountPoint *)xmalloc((size_t)npoints * sizeof(CountPoint));
    ensemble_counting_function(p_goe, N_test, npaths, burn, dt, npoints, cp,
                               E_max);

    printf("    E      N_J(E)/N   N_rvm(E)/N    diff\n");
    printf("    ─────  ─────────  ──────────    ────────\n");
    for (int i = 0; i < npoints; i += 3) {
      printf("    %5.1f  %9.4f  %10.4f   %+8.4f\n", cp[i].energy, cp[i].N_emp,
             cp[i].N_rvm, cp[i].diff);
    }
    free(cp);
  }

  /* ── Best potential from Phase 3 ── */
  section_header("Optimized V(x) = 0.2·x²·log(1+0.001·x²)");

  {
    VParams p_opt = {0.0, 0.2, 0.001};
    CountPoint *cp =
        (CountPoint *)xmalloc((size_t)npoints * sizeof(CountPoint));
    ensemble_counting_function(p_opt, N_test, npaths, burn, dt, npoints, cp,
                               E_max);

    printf("    E      N_J(E)/N   N_rvm(E)/N    diff\n");
    printf("    ─────  ─────────  ──────────    ────────\n");
    for (int i = 0; i < npoints; i += 3) {
      printf("    %5.1f  %9.4f  %10.4f   %+8.4f\n", cp[i].energy, cp[i].N_emp,
             cp[i].N_rvm, cp[i].diff);
    }
    free(cp);
  }

  /* ── Strong log tail ── */
  section_header("Strong log V(x) = 0.5·x²·log(1+0.01·x²)");

  {
    VParams p_strong = {0.0, 0.5, 0.01};
    CountPoint *cp =
        (CountPoint *)xmalloc((size_t)npoints * sizeof(CountPoint));
    ensemble_counting_function(p_strong, N_test, npaths, burn, dt, npoints, cp,
                               E_max);

    printf("    E      N_J(E)/N   N_rvm(E)/N    diff\n");
    printf("    ─────  ─────────  ──────────    ────────\n");
    for (int i = 0; i < npoints; i += 3) {
      printf("    %5.1f  %9.4f  %10.4f   %+8.4f\n", cp[i].energy, cp[i].N_emp,
             cp[i].N_rvm, cp[i].diff);
    }
    free(cp);
  }

  printf("\n  Phase 4 complete.\n");
  printf("  Full determinant identity requires Phase 4b: coefficient-level\n");
  printf(
      "  comparison of E[det(zI-H)] to xi(1/2+iz) via LeVerrier averaging.\n");

  return 0;
}
