/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   de Boor-Golub inverse spectral reconstruction
 * @paper   yamaguchi-rh-2026.tex, §6.2
 * @theorem Lemma I
 * @proof   Eigenvalue-to-entry reconstruction
 * @step    1 — Weyl asymptotics (eigenvalue approximation)
 *
 * midpoint_reconstruct.c — universal de Boor-Golub with midpoints
 * μ_k = (γ_k + γ_{k+1})/2 — always strictly interlaces.
 * Works for all N without exception. The default reconstruction.
 */

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

static const double zeta[30] = {
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062, 37.586178,
    40.918719, 43.327073, 48.005151, 49.773832, 52.970321, 56.446248,
    59.347044, 60.831779, 65.112544, 67.079811, 69.546402, 72.067158,
    75.704691, 77.144840, 79.337375, 82.910381, 84.735493, 87.425275,
    88.809111, 92.491899, 94.651344, 95.870634, 98.831194, 101.317851};

static int deboor(const double *lam, const double *mu, int N, double *a,
                  double *b) {
  for (int k = 0; k < N - 1; k++)
    if (!(lam[k] < mu[k] && mu[k] < lam[k + 1]))
      return -1;
  double w[30], ws = 0;
  for (int k = 0; k < N; k++) {
    double n = 1;
    for (int j = 0; j < N - 1; j++)
      n *= lam[k] - mu[j];
    double d = 1;
    for (int j = 0; j < N; j++)
      if (j != k)
        d *= lam[k] - lam[j];
    w[k] = n / d;
    if (w[k] < 0)
      return -2;
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
  return 0;
}

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

  // Midpoint gauge: always strictly interlaces. Paper Section 4.1, Lemma I.
  printf("Midpoint reconstruction: mu_k = (gamma_k + gamma_{k+1})/2\n");
  printf("Always interlaces — universal for all N.\n\n");
  printf(" N   fwd        a_0        a_1        a_2\n");
  printf("---  ---------  ---------  ---------  ---------\n");
  for (int N = 3; N <= 30; N++) {
    double lam[30], mu[29], a[30], b[29];
    for (int k = 0; k < N; k++)
      lam[k] = zeta[k];
    for (int k = 0; k < N - 1; k++)
      mu[k] = 0.5 * (zeta[k] + zeta[k + 1]);
    int ret = deboor(lam, mu, N, a, b);
    if (ret != 0) {
      printf("%2d  FAILED (%d)\n", N, ret);
      break;
    }

    double M[900];
    int n = N;
    for (int i = 0; i < n * n; i++)
      M[i] = 0;
    for (int k = 0; k < n; k++) {
      M[k * n + k] = a[k];
      if (k < n - 1) {
        M[k * n + k + 1] = b[k];
        M[(k + 1) * n + k] = b[k];
      }
    }
    for (int sw = 0; sw < 50; sw++) {
      double moff = 0;
      for (int p = 0; p < n - 1; p++)
        for (int q = p + 1; q < n; q++) {
          double v = fabs(M[p * n + q]);
          if (v > moff)
            moff = v;
        }
      if (moff < 1e-14)
        break;
      for (int p = 0; p < n - 1; p++)
        for (int q = p + 1; q < n; q++) {
          double apq = M[p * n + q];
          if (fabs(apq) < 1e-16 * (fabs(M[p * n + p]) + fabs(M[q * n + q]) + 1))
            continue;
          double app = M[p * n + p], aqq = M[q * n + q];
          double tau = (aqq - app) / (2 * apq);
          double t = tau >= 0 ? 1 / (tau + sqrt(1 + tau * tau))
                              : -1 / (-tau + sqrt(1 + tau * tau));
          double c = 1 / sqrt(1 + t * t), s = t * c;
          for (int i = 0; i < n; i++) {
            double vip = M[i * n + p], viq = M[i * n + q];
            M[i * n + p] = vip * c - viq * s;
            M[i * n + q] = vip * s + viq * c;
          }
          for (int j = 0; j < n; j++) {
            double vpj = M[p * n + j], vqj = M[q * n + j];
            M[p * n + j] = vpj * c - vqj * s;
            M[q * n + j] = vpj * s + vqj * c;
          }
        }
    }
    double ev[30];
    for (int k = 0; k < n; k++)
      ev[k] = M[k * n + k];
    for (int i = 1; i < n; i++) {
      double key = ev[i];
      int j = i - 1;
      while (j >= 0 && ev[j] > key) {
        ev[j + 1] = ev[j];
        j--;
      }
      ev[j + 1] = key;
    }
    double max_err = 0;
    for (int k = 0; k < N; k++)
      if (fabs(ev[k] - lam[k]) > max_err)
        max_err = fabs(ev[k] - lam[k]);
    printf("%2d  %.1e  %9.3f  %9.3f  %9.3f\n", N, max_err, a[0], a[1], a[2]);
  }
  printf("\nMidpoints: universal, always interlaces. No exceptions.\n");
  return 0;
}
