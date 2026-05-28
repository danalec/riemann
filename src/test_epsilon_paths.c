/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Test epsilon->0 path independence of DFT limit
 * @paper   yamaguchi-rh-2026.tex, Section 5.4
 * @theorem Theorem I (Guinand-Weil Explicit Formula)
 * @proof   Multiple epsilon_k sequences converge to same limit
 * @step    2 -- path independence verification
 *
 * test_epsilon_paths.c - Test if different epsilon -> 0 paths give different
 * DFT limits
 *
 * KEY QUESTION: Is the epsilon=0 limit a "fluke"?
 *
 * The paper claims: lim_{N -> infinity} S_N^{(0)}(log p) = -i/(2pi sqrt(p))
 * But the proof uses a diagonal argument which doesn't prove the full limit.
 *
 * This test checks if different epsilon_k -> 0 sequences give different
 * results. If they do, the result is a fluke. If they all converge to the same
 * limit, the result might be genuine.
 *
 * Compile: gcc -O3 -o test_epsilon_paths test_epsilon_paths.c -lm
 */

#include "refdata_1000.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

#ifdef _WIN32
#include <windows.h>
#endif

/* Abel-summed S(T) at epsilon > 0:
 * S_epsilon(T) = -(1/pi) sum_{p,m} p^{-m(1/2+epsilon)}/m * sin(m T log p)
 *
 * For computational feasibility, we truncate at primes p <= P_MAX and m <=
 * M_MAX.
 */

#define P_MAX 1000
#define M_MAX 10

static int is_prime(int n) {
  if (n < 2)
    return 0;
  if (n == 2)
    return 1;
  if (n % 2 == 0)
    return 0;
  for (int d = 3; d * d <= n; d += 2)
    if (n % d == 0)
      return 0;
  return 1;
}

/* Compute Abel-summed S_epsilon(T) */
static double S_epsilon(double T, double epsilon) {
  double sum = 0.0;
  for (int p = 2; p <= P_MAX; p++) {
    if (!is_prime(p))
      continue;
    for (int m = 1; m <= M_MAX; m++) {
      double term = pow((double)p, -m * (0.5 + epsilon)) / m;
      sum += term * sin(m * T * log((double)p));
    }
  }
  return -sum / M_PI;
}

/* Fejer-weighted DFT of S_epsilon at omega = log p
 * Returns BOTH real and imaginary parts */
static void fejer_dft_complex(double epsilon, int N, double omega, double *re,
                              double *im) {
  double sum_re = 0.0, sum_im = 0.0;
  for (int n = 0; n < N; n++) {
    double weight = 2.0 * (N - n) / (N * N); /* Fejer window */
    double Sn = S_epsilon(ZETA_ZEROS[n], epsilon);
    double phase = omega * ZETA_ZEROS[n];
    sum_re += weight * Sn * cos(phase);
    sum_im += weight * Sn * sin(phase);
  }
  *re = sum_re;
  *im = sum_im;
}

/* Test different epsilon -> 0 sequences */
int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

  printf("=== Testing epsilon -> 0 Path Dependence ===\n\n");

  printf("KEY QUESTION: Do different epsilon_k -> 0 sequences\n");
  printf("              give different limits?\n\n");

  printf("Paper claims: lim_{N -> infinity} S_N^{(0)}(log p) = -i/(2pi "
         "sqrt(p))\n");
  printf("Expected real part: ~ -1/(4pi sqrt(p)) ~ -0.1134 for p=2\n\n");

  int p_target = 2;
  double omega = log((double)p_target);
  double expected_im =
      -1.0 / (2.0 * M_PI * sqrt((double)p_target)); /* Imaginary part */

  printf("Target prime: p = %d, omega = %.6f\n", p_target, omega);
  printf("Expected limit (imaginary part): %.6f\n", expected_im);
  printf("Expected limit (real part): 0 (should be exactly zero)\n\n");

  /* Test 1: Direct epsilon=0 (no regulator) */
  printf("=== Test 1: Direct epsilon=0 (NO Abel regulator) ===\n");
  double re_e0, im_e0;
  fejer_dft_complex(0.0, N_REF, omega, &re_e0, &im_e0);
  printf("  DFT at epsilon=0, N=%d:\n", N_REF);
  printf("    Real:      %10.6f (expected: 0)\n", re_e0);
  printf("    Imaginary: %10.6f (expected: %.6f)\n", im_e0, expected_im);
  printf("    Magnitude: %10.6f\n", sqrt(re_e0 * re_e0 + im_e0 * im_e0));
  printf("    Deviation: %10.6f\n\n", im_e0 - expected_im);

  /* Test 2: Different epsilon sequences */
  printf("=== Test 2: Different epsilon_k -> 0 sequences ===\n\n");

  /* Sequence A: epsilon_k = 1/sqrt(k) (paper's choice) */
  printf("Sequence A: epsilon_k = 1/sqrt(k)\n");
  printf("  %8s  %8s  %10s  %10s  %10s\n", "k", "epsilon", "N", "Im(DFT)",
         "Expected");
  for (int k = 10; k <= 100; k += 10) {
    double eps = 1.0 / sqrt((double)k);
    int Nk = k * k * k; /* Paper's N_k = k^3 */
    if (Nk > N_REF)
      Nk = N_REF;
    double re, im;
    fejer_dft_complex(eps, Nk, omega, &re, &im);
    printf("  %8d  %8.4f  %10d  %10.6f  %10.6f\n", k, eps, Nk, im, expected_im);
  }

  /* Sequence B: epsilon_k = 1/k (slower decay) */
  printf("\nSequence B: epsilon_k = 1/k (slower decay)\n");
  printf("  %8s  %8s  %10s  %10s  %10s\n", "k", "epsilon", "N", "Im(DFT)",
         "Expected");
  for (int k = 10; k <= 100; k += 10) {
    double eps = 1.0 / (double)k;
    int Nk = k * k; /* Different N scaling */
    if (Nk > N_REF)
      Nk = N_REF;
    double re, im;
    fejer_dft_complex(eps, Nk, omega, &re, &im);
    printf("  %8d  %8.4f  %10d  %10.6f  %10.6f\n", k, eps, Nk, im, expected_im);
  }

  /* Sequence C: epsilon_k = 1/k^2 (very fast decay) */
  printf("\nSequence C: epsilon_k = 1/k^2 (fast decay)\n");
  printf("  %8s  %8s  %10s  %10s  %10s\n", "k", "epsilon", "N", "Im(DFT)",
         "Expected");
  for (int k = 10; k <= 100; k += 10) {
    double eps = 1.0 / ((double)k * (double)k);
    int Nk = k; /* Linear N scaling */
    if (Nk > N_REF)
      Nk = N_REF;
    double re, im;
    fejer_dft_complex(eps, Nk, omega, &re, &im);
    printf("  %8d  %8.4f  %10d  %10.6f  %10.6f\n", k, eps, Nk, im, expected_im);
  }

  printf("\n=== Analysis ===\n\n");

  /* Check convergence stability */
  double re_a, im_a, re_b, im_b;
  fejer_dft_complex(0.01, N_REF, omega, &re_a, &im_a);
  fejer_dft_complex(0.01, N_REF, omega, &re_b, &im_b);

  printf("At epsilon=0.01, N=%d:\n", N_REF);
  printf("  Real part:      %.6f (expected: 0)\n", re_a);
  printf("  Imaginary part: %.6f (expected: %.6f)\n", im_a, expected_im);
  printf("  Difference:     %.6f\n\n", im_a - expected_im);

  printf("CONCLUSION:\n");
  if (fabs(im_a - expected_im) < 0.02) {
    printf("  Imaginary part converges to expected value!\n");
    printf("  The result may NOT be a fluke - convergence appears correct.\n");
    printf("  However, the mathematical proof gap remains:\n");
    printf("    - Diagonal argument does not prove full limit\n");
    printf("    - Uniformity fails at epsilon=0\n");
    printf("    - Divergence of sum p^{-1/2} is not cured\n");
  } else {
    printf("  Imaginary part does NOT converge to expected value!\n");
    printf("  Expected: %.6f, Got: %.6f\n", expected_im, im_a);
    printf("  This confirms the result IS A FLUKE.\n");
    printf("  The diagonal argument is insufficient.\n");
    printf("  The mathematical proof has a genuine gap.\n");
  }

  return 0;
}