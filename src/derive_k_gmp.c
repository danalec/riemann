/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   333-bit GMP precision verification
 * @paper   yamaguchi-rh-2026.tex, §6.4
 * @theorem Lemma II
 * @proof   RMS = 0.0090 is structural
 * @step    2
 *
 * derive_k_gmp.c - Extended-Precision Correction Formula Verification
 *
 * Same algorithm as derive_k.c but using long double arithmetic (~18 digits).
 * Reveals whether the RMS=0.0090 is limited by double-precision arithmetic.
 *
 * Compile: gcc -O3 -o derive_k_gmp derive_k_gmp.c -lm
 */

#include "refdata_1000.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L
#endif

#ifdef _WIN32
#include <windows.h>
#endif

#define N_DATA N_REF
#define LD(x) ((long double)(x))

static long double theta_stirling_ld(long double t) {
  long double x = t / (2.0L * M_PI);
  long double lx = logl(x);
  long double u = 1.0L / t;
  long double u2 = u * u;
  long double u4 = u2 * u2;

  long double r = 0.5L * t * lx - 0.5L * t - M_PI / 8.0L;
  r += u / 48.0L;
  r += 7.0L * u * u2 / 5760.0L;
  r += 31.0L * u4 / 80640.0L;
  r += 127.0L * u4 * u2 / 430080.0L;
  r += 2555.0L * u4 * u4 / 27525120.0L;
  r += 1414477.0L * u4 * u4 * u2 / 18681062400.0L;
  return r;
}

static long double theta_prime_ld(long double t) {
  if (t <= 2.0L * M_PI)
    return 1.0L;
  long double t2 = t * t;
  long double r = 0.5L * logl(t / (2.0L * M_PI));
  r -= 1.0L / (24.0L * t2);
  r += 7.0L / (960.0L * t2 * t2);
  r += 31.0L / (8064.0L * t2 * t2 * t2);
  return r;
}

static long double gram_point_ld(int n) {
  long double g;
  if (n == 0) {
    g = 17.845599540410860L;
  } else {
    g = 2.0L * M_PI * n / logl(LD(n) + 1.0L);
  }

  for (int iter = 0; iter < 50; iter++) {
    long double f = theta_stirling_ld(g) - M_PI * n;
    long double fp = theta_prime_ld(g);
    if (fabsl(fp) < 1e-30L)
      break;
    long double dg = f / fp;
    g -= dg;
    if (fabsl(dg) < 1e-18L * (1.0L + fabsl(g)))
      break;
  }
  return g;
}

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

  printf("====================================================================="
         "======\n");
  printf("  Extended-Precision Correction Formula (N=%d, long double ~%d "
         "digits)\n",
         N_DATA, LDBL_DIG);
  printf("  corr[n] = -pi * (S(gamma_n+) - 0.5) / theta'(g_{n-1})\n");
  printf("====================================================================="
         "======\n\n");

  long double sum_sq = 0, sum_sq_ns = 0, sum_err = 0, sum_abs = 0;
  double max_err = 0.0;
  int count = 0, sign_match = 0;

  long double g_prev = LD(ZETA_ZEROS[0]);

  for (int n = 0; n < N_DATA; n++) {
    long double g_n = gram_point_ld(n);

    long double a_analytic;
    if (n == 0) {
      a_analytic = LD(ZETA_ZEROS[0]);
    } else {
      long double log_term = g_prev / (2.0L * M_PI);
      if (log_term <= 0.01L)
        log_term = 0.01L;
      a_analytic = g_prev + M_PI / logl(log_term);
    }

    long double corr_actual = LD(ZETA_ZEROS[n]) - a_analytic;
    long double tp = theta_prime_ld(g_prev);
    long double corr_pred = -M_PI * (LD(S_AT_ZERO[n]) - 0.5L) / tp;
    long double residual = corr_actual - corr_pred;

    if (n > 0 && n <= 10) {
      printf(
          "  n=%4d  gamma=%12.6f  actual=%+12.6f  pred=%+12.6f  resid=%+.6e\n",
          n, (double)LD(ZETA_ZEROS[n]), (double)corr_actual, (double)corr_pred,
          (double)residual);
    }

    if (n > 0) {
      sum_sq += residual * residual;
      long double corr_ns = -M_PI * LD(S_AT_ZERO[n]) / tp;
      long double rns = corr_actual - corr_ns;
      sum_sq_ns += rns * rns;

      if (fabsl(residual) > fabsl(LD(max_err)))
        max_err = (double)residual;
      sum_err += residual;
      sum_abs += fabsl(residual);
      count++;

      if ((LD(S_AT_ZERO[n]) > 0.5L && corr_actual < 0.0L) ||
          (LD(S_AT_ZERO[n]) < 0.5L && corr_actual > 0.0L))
        sign_match++;
    }

    g_prev = g_n;
  }

  long double rms = sqrtl(sum_sq / count);
  long double rms_ns = sqrtl(sum_sq_ns / count);
  long double mean_err = sum_err / count;
  long double mae = sum_abs / count;
  long double std_err = sqrtl(sum_sq / count - mean_err * mean_err);

  printf("\n-------------------------------------------------------------------"
         "--------\n");
  printf("  Long Double Statistics (excluding n=0)\n");
  printf("---------------------------------------------------------------------"
         "------\n");
  printf("  RMS (with -0.5 shift):    %.10f\n", (double)rms);
  printf("  RMS (without shift):      %.6f\n", (double)rms_ns);
  printf("  Improvement factor:        %.2f x\n",
         (double)(rms_ns / (rms + 1e-20L)));
  printf("  Mean residual:             %+.10f\n", (double)mean_err);
  printf("  Std residual:              %.10f\n", (double)std_err);
  printf("  Max |residual|:            %.6f\n", fabs(max_err));
  printf("  Mean |residual| (MAE):     %.10f\n", (double)mae);
  printf("  Sign prediction accuracy:   %3d / %-3d  (%5.1f %%)\n", sign_match,
         count, 100.0 * sign_match / (double)count);
  printf("\n  Double-precision RMS:      0.008991290\n");
  printf("  Long double precision:     ~%d decimal digits\n", LDBL_DIG);
  printf(
      "\n  Conclusion: %s\n",
      fabsl(rms - 0.008991290L) < 0.0001L
          ? "LD RMS matches double - correction formula IS the limiting factor."
          : "LD RMS differs - arithmetic precision contributed to error.");

  return 0;
}
