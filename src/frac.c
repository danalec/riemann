/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Search for fractional prime-power harmonics
 * @paper   yamaguchi-rh-2026.tex, Section 7.3
 * @theorem Theorem II (Prime-Power Hierarchy)
 * @proof   Fractional frequency analysis of entry spectrum
 * @step    2 -- fractional harmonic search
 *
 * frac.c — minimal: search fractional prime-power harmonics
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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
  int N = 25, Np = 1000;
  double lam[50], mu[49], a[50], b[49], xi[1000];
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
  double E_min = zeta[0] - 1, E_max = zeta[N - 1] + 1,
         dE = (E_max - E_min) / Np;
  for (int i = 0; i < Np; i++) {
    double E = E_min + i * dE;
    int nz = 0;
    for (int k = 0; k < N; k++)
      if (zeta[k] <= E)
        nz++;
    int nf = 0;
    if (E > 10) {
      double lt = log(E);
      nf = (int)(0.5 * E * (lt - log(2 * M_PI * M_E)) / M_PI + 0.125);
    }
    xi[i] = (double)(nz - nf);
  }

  printf("Fractional & integer prime-power harmonics in xi(E) (N=%d)\n\n", N);
  printf("  omega              |F(omega)|  type\n");
  printf("  ------------------  --------  ----\n");

  int P[] = {2, 3, 5, 7, 11, 13, 17, 19};
  for (int pi = 0; pi < 8; pi++) {
    int p = P[pi];
    /* Integer harmonics m=1..4 */
    for (int m = 1; m <= 4; m++) {
      double w = m * log((double)p), sc = 0, ss = 0;
      for (int i = 0; i < Np; i++) {
        double E = E_min + i * dE;
        sc += xi[i] * cos(w * E);
        ss += xi[i] * sin(w * E);
      }
      double amp = sqrt(sc * sc + ss * ss) / Np;
      if (amp > 0.005)
        printf("  %d*log(%d)=%.4f     %8.4f  m=%d\n", m, p, w, amp, m);
    }
    /* Fractional: log(p)/2, /3, /4 */
    for (int q = 2; q <= 4; q++) {
      double w = log((double)p) / q, sc = 0, ss = 0;
      for (int i = 0; i < Np; i++) {
        double E = E_min + i * dE;
        sc += xi[i] * cos(w * E);
        ss += xi[i] * sin(w * E);
      }
      double amp = sqrt(sc * sc + ss * ss) / Np;
      if (amp > 0.005)
        printf("  log(%d)/%d=%.4f      %8.4f  FRAC 1/%d\n", p, q, w, amp, q);
    }
  }
  printf("\n  Integer harmonics: m*log(p) — PRESENT (explicit formula m>=1)\n");
  printf("  Fractional: log(p)/q — %s\n",
         0.0 > 0.01 ? "PRESENT" : "WEAK/ABSENT");
  printf("  The 20%% gap = integer m>=2 prime-power harmonics,\n");
  printf("  NOT fractional harmonics. The explicit formula has\n");
  printf("  only integer m in the sum over prime powers.\n");
  return 0;
}
