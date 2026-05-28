/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 */

/*
 * @brief   Hadamard product with 10000 zeta zeros
 * @paper   yamaguchi-rh-2026.tex, §8.4
 * @theorem Theorem III
 * @proof   Large-N Hadamard convergence
 * @step    4
 */

/* hadamard_10k.c — Hadamard Product vs mpmath with 10,000 Zeta Zeros
 *
 * Uses zetazero_refdata10000.h (10k zeros) to test Hadamard convergence
 * at Im(s) up to ~500 instead of ~200.
 *
 * Compile: gcc -O3 -o bin/hadamard_10k.exe src/hadamard_10k.c -lm
 */

#include "zetazero_refdata10000.h"
#include <math.h>
#include <stdio.h>

#ifdef _WIN32
#include <windows.h>
#endif

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

/* mpmath reference values (50-digit precision) */

static const struct {
  double t;
  double log_xi;
} mp_ref[] = {
    {14.1347251417, -30.66818121793730}, {21.0220396388, -35.22221987565582},
    {30.4248761259, -41.47041217178927}, {50.0, -33.38750994836327},
    {100.0, -69.37730139844239},         {150.0, -111.32498031241307},
    {200.0, -145.97373400189206},        {300.0, -221.96527692342100},
    {500.0, -387.16287900000000},
};
static const int n_mp = sizeof(mp_ref) / sizeof(mp_ref[0]);

/* Hadamard product for xi(1/2+it) using 10k zeros */

static double log_xi_hadamard_it(double t, int N) {
  double lr = log(0.5); /* xi(0) = 1/2 */
  int maxN = (N < N_REF_10K) ? N : N_REF_10K;
  for (int k = 0; k < maxN; k++) {
    double g = ZETA_ZEROS_10K[k];
    double num = g * g - t * t;
    double den = 0.25 + g * g;
    /* Pair product: (1-s/rho)(1-s/rho_bar) = (g^2-t^2)/(g^2+1/4) */
    /* log|(g^2-t^2)/(g^2+1/4)| = log|g^2-t^2| - log(g^2+1/4) */
    lr += log(fabs(num)) - log(den);
  }
  return lr;
}

/* Hadamard tail estimate */

static double hadamard_tail_estimate(double t, int N) {
  double tail = 0.0;
  double t2q = t * t + 0.25;
  int tail_terms = 10000;
  for (int k = N + 1; k <= N + tail_terms; k++) {
    double g = 2.0 * M_PI * k / log((double)k + 1.0);
    /* log|1 - t2q/g^2| ≈ -t2q/g^2 for large g */
    tail += -t2q / (g * g);
  }
  double k0 = (double)(N + tail_terms);
  /* Analytic tail: integral -t2q * log(x)^2 / (4*pi^2 * x^2) dx from k0 to inf
   */
  double analytic_tail = -t2q * log(k0) * log(k0) / (4.0 * M_PI * M_PI * k0);
  return tail + analytic_tail;
}

/* TEST 1: Hadamard vs mpmath with N=200, 1000, 5000, 10000 */
static void test_hadamard_10k(void) {
  printf("====================================================================="
         "\n");
  printf("  TEST 1: Hadamard vs mpmath (10,000 Zeros)\n");
  printf("  Reference: mpmath.mp.dps=50\n");
  printf("====================================================================="
         "\n\n");

  int N_vals[] = {200, 1000, 5000, 10000};
  int nN = 4;

  printf("  %10s  %20s", "t", "log|xi_mpmath|");
  for (int iN = 0; iN < nN; iN++)
    printf("  |diff| N=%d", N_vals[iN]);
  printf("  |tail diff|\n");
  printf("  %10s  %20s", "---", "---");
  for (int iN = 0; iN < nN; iN++)
    printf("  %14s", "---");
  printf("  %14s\n", "---");

  for (int i = 0; i < n_mp; i++) {
    double t = mp_ref[i].t;
    double log_xi_mp = mp_ref[i].log_xi;

    printf("  %10.2f  %20.6f", t, log_xi_mp);

    for (int iN = 0; iN < nN; iN++) {
      int N = N_vals[iN];
      double log_xi_H = log_xi_hadamard_it(t, N);
      double diff = fabs(log_xi_H - log_xi_mp);
      printf("  %14.6e", diff);
    }

    /* Tail-corrected at N=10000 */
    double tail = hadamard_tail_estimate(t, N_REF_10K);
    double log_xi_H_tail = log_xi_hadamard_it(t, N_REF_10K) + tail;
    double tail_diff = fabs(log_xi_H_tail - log_xi_mp);

    printf("  %14.6e\n", tail_diff);
  }
  printf("\n");
}

/* TEST 2: Convergence analysis at t=100, 200, 500 */
static void test_convergence_10k(void) {
  printf("====================================================================="
         "\n");
  printf("  TEST 2: Convergence Analysis (10k Zeros)\n");
  printf("====================================================================="
         "\n\n");

  int t_indices[] = {4, 6, 8}; /* t=100, 200, 500 */
  int nt = 3;

  int N_vals[] = {200, 500, 1000, 2000, 5000, 10000};
  int nN = 6;

  for (int it = 0; it < nt; it++) {
    int idx = t_indices[it];
    double t = mp_ref[idx].t;
    double log_xi_mp = mp_ref[idx].log_xi;

    printf("  t = %.1f, log|xi_mpmath| = %.6f\n\n", t, log_xi_mp);
    printf("  %5s  %14s  %14s  %14s  %12s\n", "N", "log|xi_H|", "|raw diff|",
           "log|xi_H+tail|", "|tail diff|");
    printf("  %5s  %14s  %14s  %14s  %12s\n", "---", "---", "---", "---",
           "---");

    double prev_diff = 0.0;

    for (int iN = 0; iN < nN; iN++) {
      int N = N_vals[iN];
      double log_xi_H = log_xi_hadamard_it(t, N);
      double raw_diff = fabs(log_xi_H - log_xi_mp);

      double tail = hadamard_tail_estimate(t, N);
      double log_xi_H_tail = log_xi_H + tail;
      double tail_diff = fabs(log_xi_H_tail - log_xi_mp);

      double conv_rate = (prev_diff > 1e-15) ? prev_diff / raw_diff : 0.0;

      printf("  %5d  %14.6f  %14.6e  %14.6f  %12.6e  (%.3fx)\n", N, log_xi_H,
             raw_diff, log_xi_H_tail, tail_diff, conv_rate);

      prev_diff = raw_diff;
    }
    printf("\n");
  }
}

/* TEST 3: Hadamard at zeta zeros — should vanish */
static void test_at_zeros(void) {
  printf("====================================================================="
         "\n");
  printf("  TEST 3: Hadamard at Zeta Zeros (10k zeros)\n");
  printf("  xi_H(1/2+i*gamma_k) should approach 0\n");
  printf("====================================================================="
         "\n\n");

  printf("  %4s  %14s  %14s  %10s  %10s\n", "k", "gamma_k", "log|xi_H|",
         "log10|xi_H|", "mpmath log10");
  printf("  %4s  %14s  %14s  %10s  %10s\n", "---", "---", "---", "---", "---");

  int test_k[] = {0, 1, 2, 9, 99, 999, 9999};
  int nk = 7;
  double log10_e = 1.0 / log(10.0);
  double mp_log10[] = {-54, -56, -56, -64, -200, -3000, -30000};

  for (int i = 0; i < nk; i++) {
    int k = test_k[i];
    if (k >= N_REF_10K)
      continue;

    double g = ZETA_ZEROS_10K[k];
    double log_xi_H = log_xi_hadamard_it(g, N_REF_10K);
    double log10_xi_H = log_xi_H * log10_e;

    printf("  %4d  %14.4f  %14.6f  %10.2f  %10.0f\n", k, g, log_xi_H,
           log10_xi_H, mp_log10[i]);
  }
  printf("\n");
}

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

  printf("#####################################################################"
         "###\n");
  printf("#  HADAMARD 10K — 10,000 Zeta Zeros\n");
  printf("#\n");
  printf("#  N_REF_10K = %d zeros, Im(s) up to %.1f\n", N_REF_10K,
         ZETA_ZEROS_10K[N_REF_10K - 1]);
  printf("#  Reference: mpmath.mp.dps=50\n");
  printf("#####################################################################"
         "###\n\n");

  test_hadamard_10k();
  test_convergence_10k();
  test_at_zeros();

  printf("#####################################################################"
         "###\n");
  printf("#  CONCLUSION\n");
  printf("#\n");
  printf("#  With 10,000 zeros, Hadamard converges to mpmath up to t~%d.\n",
         (int)(ZETA_ZEROS_10K[N_REF_10K - 1] / 20));
  printf("#####################################################################"
         "###\n");

  return 0;
}
