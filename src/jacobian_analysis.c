/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Jacobian of eigenvalue map
 * @paper   yamaguchi-rh-2026.tex, Section 7.5
 * @theorem Lemma I
 * @proof   Entry vs spectral extraction
 * @step    1 — Weyl asymptotics (eigenvalue approximation)
 *
 * jacobian_analysis.c â€” is the prime perturbation direction reachable?
 * Computes the Jacobian âˆ‚{a,b}/âˆ‚Î¼ and projects Î´a_prime, Î´b_prime
 * onto the tangent space of the isospectral manifold.
 * R² > 0.99 means the prime structure IS on the manifold.
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
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062, 37.586178, 40.918719,
    43.327073, 48.005151, 49.773832, 52.970321, 56.446248, 59.347044, 60.831779,
    65.112544, 67.079811, 69.546402, 72.067158, 75.704691, 77.144840, 79.337375,
    82.910381, 84.735493, 87.425275, 88.809111};

static int deboor(const double *lam, const double *mu, int N, double *a,
                  double *b) {
  for (int k = 0; k < N - 1; k++)
    if (!(lam[k] < mu[k] && mu[k] < lam[k + 1]))
      return -1;
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

static void prime_target(int N, double *da, double *db) {
  int pr[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
  for (int k = 0; k < N; k++) {
    da[k] = 0;
    if (k < N - 1)
      db[k] = 0;
    for (int pi = 0; pi < 10; pi++) {
      double w = log(pr[pi]);
      double a = -log(pr[pi]) / (2 * M_PI * sqrt(pr[pi]));
      da[k] += a * sin(w * k);
      if (k < N - 1)
        db[k] += a * cos(w * k);
    }
  }
}

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

  int N = 9, dim = N - 1;
  double lam[30], mu_mid[29];
  for (int k = 0; k < N; k++)
    lam[k] = zeta[k];
  for (int k = 0; k < N - 1; k++)
    mu_mid[k] = 0.5 * (zeta[k] + zeta[k + 1]);

  double a0[30], b0[29];
  deboor(lam, mu_mid, N, a0, b0);
  double da_target[50], db_target[49];
  prime_target(N, da_target, db_target);

  /* Jacobian J[17][8]: rows=entries(a0..a8,b0..b7), cols=mu0..mu7 */
  double J[50][29] = {{0}};
  double eps = 0.001;
  for (int k = 0; k < dim; k++) {
    double mu_p[29];
    memcpy(mu_p, mu_mid, (size_t)(N - 1) * sizeof(double));
    mu_p[k] += eps;
    if (!(lam[k] < mu_p[k] && mu_p[k] < lam[k + 1]))
      continue;
    double ap[30], bp[29];
    deboor(lam, mu_p, N, ap, bp);
    for (int i = 0; i < N; i++)
      J[i][k] = (ap[i] - a0[i]) / eps;
    for (int i = 0; i < N - 1; i++)
      J[N + i][k] = (bp[i] - b0[i]) / eps;
  }

  /* Target vector */
  double t[50];
  for (int i = 0; i < N; i++)
    t[i] = da_target[i];
  for (int i = 0; i < N - 1; i++)
    t[N + i] = db_target[i];

  int m = 2 * N - 1;
  double A[29][29] = {{0}}, b[29] = {0};
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int r = 0; r < m; r++)
        A[i][j] += J[r][i] * J[r][j];
    }
    for (int r = 0; r < m; r++)
      b[i] += J[r][i] * t[r];
  }
  for (int k = 0; k < dim; k++) {
    if (fabs(A[k][k]) < 1e-15) {
      for (int i = k + 1; i < dim; i++)
        if (fabs(A[i][k]) > fabs(A[k][k])) {
          for (int j = k; j < dim; j++) {
            double tmp = A[k][j];
            A[k][j] = A[i][j];
            A[i][j] = tmp;
          }
          double tmp = b[k];
          b[k] = b[i];
          b[i] = tmp;
        }
    }
    for (int i = k + 1; i < dim; i++) {
      double f = A[i][k] / A[k][k];
      for (int j = k; j < dim; j++)
        A[i][j] -= f * A[k][j];
      b[i] -= f * b[k];
    }
  }
  double x[29] = {0};
  for (int i = dim - 1; i >= 0; i--) {
    double s = b[i];
    for (int j = i + 1; j < dim; j++)
      s -= A[i][j] * x[j];
    x[i] = s / A[i][i];
  }
  double rss = 0, tss = 0;
  for (int r = 0; r < m; r++) {
    double pred = 0;
    for (int k = 0; k < dim; k++)
      pred += J[r][k] * x[k];
    double resid = t[r] - pred;
    rss += resid * resid;
    tss += t[r] * t[r];
  }

  printf("Prime direction on isospectral manifold (N=%d)\n\n", N);
  printf("RÂ² = %.6f   (%s)\n", 1.0 - rss / tss,
         1.0 - rss / tss > 0.99  ? "PRIME STRUCTURE IS ON THE MANIFOLD"
         : 1.0 - rss / tss > 0.5 ? "partially reachable"
                                 : "mostly unreachable");
  printf("Residual RMS = %.4f\n\n", sqrt(rss / m));
  printf("Unreachable component:\n  k   Î´a_resid   Î´b_resid\n");
  for (int i = 0; i < N; i++) {
    double pred = 0;
    for (int k = 0; k < dim; k++)
      pred += J[i][k] * x[k];
    printf("  %d   %+9.6f  %+9.6f\n", i, t[i] - pred,
           i < N - 1 ? t[N + i] - (0) : 0.0);
  }
  /* b residual */
  for (int i = 0; i < N - 1; i++) {
    double pred = 0;
    for (int k = 0; k < dim; k++)
      pred += J[N + i][k] * x[k];
    printf("  b%d  %9s  %+9.6f\n", i, "", t[N + i] - pred);
  }

  printf("\n100%% fit is %s achievable.\n",
         1.0 - rss / tss > 0.9999 ? "FULLY" : "NOT fully");
  return 0;
}
