/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Test gauge potential theta_k for 100% prime fit
 * @paper   yamaguchi-rh-2026.tex, Section 7.3
 * @theorem Lemma III (Isospectral Gauge Freedom)
 * @proof   Jacobian analysis of mu + theta gauge parameters
 * @step    4 -- gauge potential sufficiency test
 *
 * gauge_fit_100.c — Test if gauge potential θ_k can achieve 100% prime fit.
 * Gauge: 8 μ_k (from real DBG) + 9 θ_k (phase potentials) = 17 params.
 * Entries: 9 a_k + 8 Re(b_k) + 8 Im(b_k) = 25 real entries.
 * Compute R² via Jacobian analysis.
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
    65.112544, 67.079811, 69.546402, 72.067158, 75.704691, 77.144840};

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

static void prime_tgt(int N, double *da, double *db_re, double *db_im) {
  int pr[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
  for (int k = 0; k < N; k++) {
    da[k] = 0;
    if (k < N - 1) {
      db_re[k] = 0;
      db_im[k] = 0;
    }
    for (int pi = 0; pi < 10; pi++) {
      double w = log(pr[pi]);
      double a = -log(pr[pi]) / (2 * M_PI * sqrt(pr[pi]));
      da[k] += a * sin(w * k);
      if (k < N - 1) {
        db_re[k] += a * cos(w * k);
        db_im[k] += a * sin(w * k) * 0.5;
      }
    }
  }
}

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif
  int N = 9, dim_mu = N - 1, dim_th = N, dim_total = dim_mu + dim_th;
  int m_entries = N + 2 * (N - 1); /* a_k + Re(b_k) + Im(b_k) */
  double lam[50], mu0[49];
  for (int k = 0; k < N; k++)
    lam[k] = zeta[k];
  for (int k = 0; k < N - 1; k++)
    mu0[k] = 0.5 * (zeta[k] + zeta[k + 1]);

  double a0[50], b0[49];
  real_dbg(lam, mu0, N, a0, b0);
  double th0[50] = {0}; /* start with zero phase */

  printf("Gauge potential R² analysis (N=%d)\n", N);
  printf("Gauge: %d μ_k + %d θ_k = %d params\n", dim_mu, dim_th, dim_total);
  printf("Entries: %d a_k + %d Re(b) + %d Im(b) = %d\n\n", N, N - 1, N - 1,
         m_entries);

  /* Build Jacobian: [m_entries × dim_total] */
  double J[100][100] = {{0}};
  double eps = 1e-4;

  /* μ derivatives */
  for (int k = 0; k < dim_mu; k++) {
    double mu_p[49];
    memcpy(mu_p, mu0, (size_t)dim_mu * sizeof(double));
    mu_p[k] += eps;
    if (!(lam[k] < mu_p[k] && mu_p[k] < lam[k + 1]))
      continue;
    double ap[50], bp[49];
    real_dbg(lam, mu_p, N, ap, bp);
    for (int i = 0; i < N; i++)
      J[i][k] = (ap[i] - a0[i]) / eps;
    for (int i = 0; i < N - 1; i++) {
      J[N + i][k] = (bp[i] - b0[i]) / eps;
      J[N + (N - 1) + i][k] = 0;
    }
  }

  /* θ derivatives */
  for (int k = 0; k < dim_th; k++) {
    /* Perturb θ_k: changes Re(b_k) and Im(b_k), a_k unchanged */
    double delta = eps;
    /* Effect on b_{k-1}: phase diff θ_{k-1}-θ_k changes */
    if (k > 0 && k - 1 < N - 1) {
      double ph0 = th0[k - 1] - th0[k], ph1 = ph0 - delta;
      J[N + (k - 1)][dim_mu + k] = (b0[k - 1] * (cos(ph1) - cos(ph0))) / eps;
      J[N + (N - 1) + (k - 1)][dim_mu + k] =
          (b0[k - 1] * (sin(ph1) - sin(ph0))) / eps;
    }
    /* Effect on b_k: phase diff θ_k-θ_{k+1} changes */
    if (k < N - 1) {
      double ph0 = th0[k] - th0[k + 1], ph1 = ph0 + delta;
      J[N + k][dim_mu + k] = (b0[k] * (cos(ph1) - cos(ph0))) / eps;
      J[N + (N - 1) + k][dim_mu + k] = (b0[k] * (sin(ph1) - sin(ph0))) / eps;
    }
  }

  /* Target vector */
  double t[100], da[50], dbr[50], dbi[50];
  prime_tgt(N, da, dbr, dbi);
  for (int i = 0; i < N; i++)
    t[i] = da[i];
  for (int i = 0; i < N - 1; i++) {
    t[N + i] = dbr[i];
    t[N + (N - 1) + i] = dbi[i];
  }

  /* Normal equations J^T·J·x = J^T·t */
  double A[100][100] = {{0}}, b[100] = {0};
  for (int i = 0; i < dim_total; i++) {
    for (int j = 0; j < dim_total; j++) {
      for (int r = 0; r < m_entries; r++)
        A[i][j] += J[r][i] * J[r][j];
    }
    for (int r = 0; r < m_entries; r++)
      b[i] += J[r][i] * t[r];
  }
  for (int k = 0; k < dim_total; k++) {
    if (fabs(A[k][k]) < 1e-15)
      for (int i = k + 1; i < dim_total; i++)
        if (fabs(A[i][k]) > fabs(A[k][k])) {
          for (int j = k; j < dim_total; j++) {
            double tmp = A[k][j];
            A[k][j] = A[i][j];
            A[i][j] = tmp;
          }
          double tmp = b[k];
          b[k] = b[i];
          b[i] = tmp;
        }
    for (int i = k + 1; i < dim_total; i++) {
      double f = A[i][k] / A[k][k];
      for (int j = k; j < dim_total; j++)
        A[i][j] -= f * A[k][j];
      b[i] -= f * b[k];
    }
  }
  double x[100] = {0};
  for (int i = dim_total - 1; i >= 0; i--) {
    double s = b[i];
    for (int j = i + 1; j < dim_total; j++)
      s -= A[i][j] * x[j];
    x[i] = s / A[i][i];
  }
  double rss = 0, tss = 0;
  for (int r = 0; r < m_entries; r++) {
    double pred = 0;
    for (int k = 0; k < dim_total; k++)
      pred += J[r][k] * x[k];
    double resid = t[r] - pred;
    rss += resid * resid;
    tss += t[r] * t[r];
  }

  double R2 = 1.0 - rss / tss;
  printf("R² = %.6f (%s)\n", R2,
         R2 > 0.9999 ? "100%% ACHIEVED!"
         : R2 > 0.95 ? "NEARLY"
         : R2 > 0.8  ? "partial"
                     : "low");
  printf("Rank of Jacobian: %d params for %d entries\n", dim_total, m_entries);

  /* Breakdown by channel */
  double rss_a = 0, rss_bre = 0, rss_bim = 0, tss_a = 0, tss_bre = 0,
         tss_bim = 0;
  for (int i = 0; i < N; i++) {
    double p = 0;
    for (int k = 0; k < dim_total; k++)
      p += J[i][k] * x[k];
    rss_a += (t[i] - p) * (t[i] - p);
    tss_a += t[i] * t[i];
  }
  for (int i = 0; i < N - 1; i++) {
    double p = 0;
    for (int k = 0; k < dim_total; k++)
      p += J[N + i][k] * x[k];
    rss_bre += (t[N + i] - p) * (t[N + i] - p);
    tss_bre += t[N + i] * t[N + i];
  }
  for (int i = 0; i < N - 1; i++) {
    double p = 0;
    for (int k = 0; k < dim_total; k++)
      p += J[N + N - 1 + i][k] * x[k];
    rss_bim += (t[N + N - 1 + i] - p) * (t[N + N - 1 + i] - p);
    tss_bim += t[N + N - 1 + i] * t[N + N - 1 + i];
  }
  printf("  a_k channel R²: %.4f\n", tss_a > 0 ? 1 - rss_a / tss_a : 0);
  printf("  Re(b) channel R²: %.4f\n", tss_bre > 0 ? 1 - rss_bre / tss_bre : 0);
  printf("  Im(b) channel R²: %.4f\n", tss_bim > 0 ? 1 - rss_bim / tss_bim : 0);
  return 0;
}
