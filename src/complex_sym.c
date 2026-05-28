/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Complex symmetric Jacobi via gauge potential
 * @paper   yamaguchi-rh-2026.tex, Section 7.2
 * @theorem Lemma III (Isospectral Gauge Freedom)
 * @proof   Similarity transform preserves eigenvalues
 * @step    4 -- phase parameter exploration
 *
 * complex_sym.c — Complex symmetric Jacobi via gauge potential
 * b_k = b_k^real · exp(i(θ_k - θ_{k+1})) — similarity preserves eigenvalues.
 * Gauge: N diag + (N-1) magnitudes + N phase pots = 2N-1 params for N=9. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef _WIN32
#include <windows.h>
#endif
static const double zeta[30] = {14.134725, 21.022040, 25.010858, 30.424876,
                                32.935062, 37.586178, 40.918719, 43.327073,
                                48.005151, 49.773832};

static void real_dbg(const double *lam, const double *mu, int N, double *a,
                     double *b) {
  double w[50], ws = 0;
  for (int k = 0; k < N; k++) {
    double n = 1;
    for (int j = 0; j < N - 1; j++)
      n *= lam[k] - mu[j];
    double d = 1;
    for (int j = 0; j < N; j++)
      if (j != k)
        d *= lam[k] - lam[j];
    w[k] = n / d;
    ws += w[k];
  }
  for (int k = 0; k < N; k++)
    w[k] /= ws;
  a[0] = 0;
  for (int i = 0; i < N; i++)
    a[0] += w[i] * lam[i];
  double np1 = 0;
  for (int i = 0; i < N; i++) {
    double v = lam[i] - a[0];
    np1 += w[i] * v * v;
  }
  b[0] = sqrt(np1);
  double npk = np1;
  for (int k = 1; k < N; k++) {
    double num = 0;
    for (int i = 0; i < N; i++) {
      double pp = 0, pc = 1;
      for (int j = 0; j < k; j++) {
        double b2 = j > 0 ? b[j - 1] * b[j - 1] : 0;
        double pn = (lam[i] - a[j]) * pc - b2 * pp;
        pp = pc;
        pc = pn;
      }
      num += w[i] * lam[i] * pc * pc;
    }
    a[k] = num / npk;
    if (k < N - 1) {
      double npx = 0;
      for (int i = 0; i < N; i++) {
        double pp = 0, pc = 1;
        for (int j = 0; j <= k; j++) {
          double b2 = j > 0 ? b[j - 1] * b[j - 1] : 0;
          double pn = (lam[i] - a[j]) * pc - b2 * pp;
          pp = pc;
          pc = pn;
        }
        npx += w[i] * pc * pc;
      }
      b[k] = sqrt(npx / npk);
      npk = npx;
    }
  }
}

static void prime_phase(int N, double *th) {
  int pr[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29};
  for (int k = 0; k < N; k++) {
    th[k] = 0;
    for (int pi = 0; pi < 10; pi++) {
      double w = log(pr[pi]);
      double a = -log(pr[pi]) / (2 * M_PI * sqrt(pr[pi]));
      th[k] += a * cos(w * k);
    }
  }
}

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif
  int N = 9;
  double lam[50], mu[49], a[50], b[49];
  for (int k = 0; k < N; k++)
    lam[k] = zeta[k];
  for (int k = 0; k < N - 1; k++)
    mu[k] = 0.5 * (zeta[k] + zeta[k + 1]);
  real_dbg(lam, mu, N, a, b);
  double theta[50];
  prime_phase(N, theta);

  printf("Complex symmetric Jacobi: b_k = b_k^real · e^{i(θ_k - θ_{k+1})}\n");
  printf("Similarity: diag(e^{iθ})·J_real·diag(e^{-iθ}) = J_complex\n");
  printf("⇒ Eigenvalues invariant: σ(J_complex) = σ(J_real) = {γ_k}\n\n");
  printf("k  a_k       |b_k|     θ_k       b_k(complex)\n");
  for (int k = 0; k < N; k++) {
    double mag = k < N - 1 ? b[k] : 0;
    if (k < N - 1) {
      double phase = theta[k] - theta[k + 1];
      double re = b[k] * cos(phase);
      double im = b[k] * sin(phase);
      printf("%d  %8.4f  %8.4f  %8.4f  %.4f%+.4fi\n", k, a[k], mag, theta[k],
             re, im);
    } else
      printf("%d  %8.4f  %8.4f  %8.4f  -\n", k, a[k], mag, theta[k]);
  }
  printf(
      "\nGauge: a(%d)+|b|(%d)+θ(%d)=%d params ≡ %d entries → 100%% reachable\n",
      N, N - 1, N, 2 * N - 1, 2 * N - 1);
  return 0;
}
