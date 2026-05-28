/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Test fractional prime-power harmonics at half-integer frequencies
 * @paper   yamaguchi-rh-2026.tex, Section 7.3
 * @theorem Theorem II (Prime-Power Hierarchy)
 * @proof   Frequency analysis at log(p)/2, log(p)/3, log(p)/4
 * @step    2 -- fractional frequency search
 *
 * fractional_test.c — Search for FRACTIONAL prime-power harmonics.
 * Tests frequencies: log(p)/2, log(p)/3, log(p)/4, ...
 * These would come from half-integer exponents in the functional equation
 * or from the gamma factor's asymptotic expansion.
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

static const double zeta[30] = {
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062, 37.586178, 40.918719,
    43.327073, 48.005151, 49.773832, 52.970321, 56.446248, 59.347044, 60.831779,
    65.112544, 67.079811, 69.546402, 72.067158, 75.704691, 77.144840, 79.337375,
    82.910381, 84.735493, 87.425275, 88.809111};

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif
  int N = 25;
  double lam[50], mu[49], a[50], b[49];
  for (int k = 0; k < N; k++)
    lam[k] = zeta[k];
  for (int k = 0; k < N - 1; k++)
    mu[k] = 0.5 * (zeta[k] + zeta[k + 1]);

  {
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

  int n_pts = 1000;
  double xi[1000], E_min = zeta[0] - 1, E_max = zeta[N - 1] + 1,
                   dE = (E_max - E_min) / (n_pts - 1);
  for (int i = 0; i < n_pts; i++) {
    double E = E_min + i * dE;
    int nz = 0;
    for (int k = 0; k < N; k++)
      if (zeta[k] <= E)
        nz++;
    int nf = 0;
    if (E > 10) {
      double lt = log(E);
      nf = (int)(0.5 * E * (lt - log(2 * M_PI * M_E)) / M_PI - 0.125);
    }
    xi[i] = (double)(nz - nf);
  }

  printf("Fractional Prime-Power Harmonics Search (N=%d)\n\n", N);
  printf("  Frequency           |F(ω)|    interpretation\n");
  printf("  -------------------  --------  --------------\n");

  int primes[] = {2, 3, 5, 7, 11, 13};
  double best_frac = 0;
  int best_p = 2, best_q = 1;

  /* Search fractional frequencies: log(p)/q for small p,q */
  for (int pi = 0; pi < 6; pi++) {
    for (int q = 2; q <= 6; q++) {
      double omega = log((double)primes[pi]) / (double)q;
      double sc = 0, ss = 0;
      for (int i = 0; i < n_pts; i++) {
        double E = E_min + i * dE;
        sc += xi[i] * cos(omega * E);
        ss += xi[i] * sin(omega * E);
      }
      double amp = sqrt(sc * sc + ss * ss) / n_pts;
      printf("  log(%d)/%d = %.4f     %8.4f  %s\n", primes[pi], q, omega, amp,
             amp > 0.03 ? "SIGNAL" : "noise");
      if (amp > best_frac) {
        best_frac = amp;
        best_p = primes[pi];
        best_q = q;
      }
    }
  }

  /* Gamma-factor related frequencies */
  printf("\n  --- Gamma-factor related frequencies ---\n");
  double gammas[] = {log(2 * M_PI) / 2.0, log(M_PI) / 2.0, log(2.0) / 2.0,
                     log(2 * M_PI) / 4.0};
  const char *names[] = {"log(2pi)/2", "log(pi)/2", "log(2)/2", "log(2pi)/4"};
  for (int gi = 0; gi < 4; gi++) {
    double omega = gammas[gi];
    double sc = 0, ss = 0;
    for (int i = 0; i < n_pts; i++) {
      double E = E_min + i * dE;
      sc += xi[i] * cos(omega * E);
      ss += xi[i] * sin(omega * E);
    }
    double amp = sqrt(sc * sc + ss * ss) / n_pts;
    printf("  %s=%.4f        %8.4f  %s\n", names[gi], omega, amp,
           amp > 0.02 ? "SIGNAL" : "noise");
  }

  printf("\n  Strongest fractional: log(%d)/%d, |F|=%.4f\n", best_p, best_q,
         best_frac);
  printf("  Fractional harmonics are %s.\n",
         best_frac > 0.03 ? "PRESENT" : "WEAK/ABSENT");
  return 0;
}
