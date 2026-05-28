/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Forward construction: primes to Jacobi matrix to eigenvalues
 * @paper   yamaguchi-rh-2026.tex, Section 5.2
 * @theorem Theorem I (Guinand-Weil Explicit Formula)
 * @proof   Prime-correlated entries produce zeta-zero eigenvalues
 * @step    1 -- forward explicit formula test
 *
 * forward_prime.c — FORWARD: primes → Jacobi matrix → eigenvalues
 * Build J from prime-correlated entries, compute eigenvalues,
 * compare to zeta zeros. Tests the explicit formula forward.
 *
 * J = J_free + J_prime where:
 *   J_free: a_k = g_k (Gram-like), b_k = 1.0 (constant coupling)
 *   J_prime: δa_k = Σ α_p·sin(log p·k), δb_k = Σ α_p·cos(log p·k)
 *
 * Hypothesis: eigenvalues of J should approach zeta zeros.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.71828182845904523536
#endif

#ifdef _WIN32
#include <windows.h>
#endif

static double gram_point(int k) {
  if (k == 0)
    return 17.84559954;
  if (k == 1)
    return 23.17028270;
  double t = 2.0 * M_PI * k / log((double)k);
  for (int it = 0; it < 8; it++) {
    double lt = log(t);
    double th = 0.5 * t * (lt - log(2.0 * M_PI * M_E)) + M_PI / 8.0;
    t -= (th - M_PI * k) / (0.5 * (lt - log(2.0 * M_PI * M_E) + 1.0));
  }
  return t;
}

static const double zeta[30] = {
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062, 37.586178,
    40.918719, 43.327073, 48.005151, 49.773832, 52.970321, 56.446248,
    59.347044, 60.831779, 65.112544, 67.079811, 69.546402, 72.067158,
    75.704691, 77.144840, 79.337375, 82.910381, 84.735493, 87.425275,
    88.809111, 92.491899, 94.651344, 95.870634, 98.831194, 101.317851};

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif
  int Nmax = 25;
  int primes[] = {2,  3,  5,  7,  11, 13, 17, 19, 23, 29, 31, 37, 41,
                  43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
  int np = 25;

  printf("FORWARD: Primes → Jacobi Matrix → Eigenvalues\n");
  printf("J = J_free + Σ α_p·P_p  with α_p = -log(p)/(2π√p)\n\n");
  printf("N   max|λ-γ|    best_match   at_k\n");
  printf("---  ----------  ----------  ----\n");

  for (int N = 3; N <= Nmax; N++) {
    double a[50], b[49], M[2500] = {0};
    /* Free Jacobi: a_k = gram_point(k), b_k = (spacing between grams)/2 */
    for (int k = 0; k < N; k++) {
      double gk = gram_point(k), gkp1 = gram_point(k + 1);
      a[k] = gk;
      if (k < N - 1)
        b[k] = 0.5 * (gkp1 - gk);
    }
    /* Add prime perturbation */
    for (int k = 0; k < N; k++) {
      double da = 0, db = 0;
      for (int pi = 0; pi < np; pi++) {
        double w = log((double)primes[pi]);
        double alpha =
            -log((double)primes[pi]) / (2.0 * M_PI * sqrt((double)primes[pi]));
        da += alpha * sin(w * k);
        if (k < N - 1)
          db += alpha * cos(w * k);
      }
      a[k] += da;
      if (k < N - 1)
        b[k] += db;
    }
    /* Build dense matrix */
    for (int k = 0; k < N; k++) {
      M[k * N + k] = a[k];
      if (k < N - 1) {
        M[k * N + k + 1] = b[k];
        M[(k + 1) * N + k] = b[k];
      }
    }
    /* Jacobi eigenvalues */
    for (int sw = 0; sw < 100; sw++) {
      double moff = 0;
      for (int p = 0; p < N - 1; p++)
        for (int q = p + 1; q < N; q++) {
          double v = fabs(M[p * N + q]);
          if (v > moff)
            moff = v;
        }
      if (moff < 1e-14)
        break;
      for (int p = 0; p < N - 1; p++)
        for (int q = p + 1; q < N; q++) {
          double apq = M[p * N + q];
          if (fabs(apq) < 1e-16 * (fabs(M[p * N + p]) + fabs(M[q * N + q]) + 1))
            continue;
          double app = M[p * N + p], aqq = M[q * N + q],
                 tau = (aqq - app) / (2 * apq),
                 t = tau >= 0 ? 1 / (tau + sqrt(1 + tau * tau))
                              : -1 / (-tau + sqrt(1 + tau * tau)),
                 c = 1 / sqrt(1 + t * t), s = t * c;
          for (int i = 0; i < N; i++) {
            double vp = M[i * N + p], vq = M[i * N + q];
            M[i * N + p] = vp * c - vq * s;
            M[i * N + q] = vp * s + vq * c;
          }
          for (int j = 0; j < N; j++) {
            double vp = M[p * N + j], vq = M[q * N + j];
            M[p * N + j] = vp * c - vq * s;
            M[q * N + j] = vp * s + vq * c;
          }
        }
    }
    double ev[50];
    for (int k = 0; k < N; k++)
      ev[k] = M[k * N + k];
    for (int i = 1; i < N; i++) {
      double key = ev[i];
      int j = i - 1;
      while (j >= 0 && ev[j] > key) {
        ev[j + 1] = ev[j];
        j--;
      }
      ev[j + 1] = key;
    }

    double max_err = 0, best_dist = 1e300;
    int best_k = 0;
    for (int k = 0; k < N; k++) {
      double d = fabs(ev[k] - zeta[k]);
      if (d > max_err)
        max_err = d;
      if (d < best_dist) {
        best_dist = d;
        best_k = k;
      }
    }
    printf("%2d  %.1e     k=%d:%.4f    k=%d\n", N, max_err, best_k, best_dist,
           best_k);
  }

  printf("\n  The forward construction (primes→matrix→eigenvalues)\n");
  printf("  produces eigenvalues that diverge from zeta zeros.\n");
  printf("  The prime-correlated perturbation alone does NOT\n");
  printf("  generate zeta zeros — it needs the spectral shift ξ(E)\n");
  printf("  to connect entries to eigenvalue positions.\n");
  printf("  This is what the de Boor-Golub inverse solves.\n");
  return 0;
}
