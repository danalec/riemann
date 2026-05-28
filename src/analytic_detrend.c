/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Analytic detrend of dBG entries
 * @paper   yamaguchi-rh-2026.tex, §6.2
 * @theorem Lemma I
 * @proof   Exploratory
 * @step    1 — Weyl asymptotics (eigenvalue approximation)
 *
 * analytic_detrend.c -- Extract prime-frequency structure from dBG entries
 *
 * PROBLEM: Raw (a_zeta - a_weyl) is O(3-7), but prime perturbation is O(0.1).
 * The dBG reconstruction is nonlinear -- small eigenvalue shifts cause large
 * entry changes. The prime signal is buried under the smooth trend.
 *
 * SOLUTION: Detrend the entries before fitting.
 * 1. Compute dBG entries for zeta zeros at N=10..50
 * 2. Fit a smooth polynomial (degree 3-5) to a_k and b_k vs k
 * 3. Compute residuals: r_a[k] = a_k - smooth_a(k)
 * 4. Fit residuals to sum_p c_p * sin(log p * k), sum_p c_p * cos(log p * k)
 * 5. Compare c_p to predicted alpha_p = -log(p)/(2pi*sqrt(p))
 * 6. Also try: detrend against Weyl-law reference, then detrend residual
 *
 * Uses Euler form: sin(k ln p) = [p^{ik} - p^{-ik}]/2i
 *                  cos(k ln p) = [p^{ik} + p^{-ik}]/2
 *                  DFT: F(omega) = (1/N) sum x_k * e^{-i*omega*k}
 *
 * Compile: gcc -Wall -Wextra -Wconversion -Wshadow -Werror -O3
 *          -fno-strict-aliasing -fno-peel-loops -fno-unswitch-loops
 *          -Isrc -o analytic_detrend src/analytic_detrend.c -lm
 */

#include "random_matrix_utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef _WIN32
#include <windows.h>
#endif

static const double ZZ[50] = {
    14.134725,  21.022040,  25.010858,  30.424876,  32.935062,  37.586178,
    40.918719,  43.327073,  48.005151,  49.773832,  52.970321,  56.446248,
    59.347044,  60.831779,  65.112544,  67.079811,  69.546402,  72.067158,
    75.704691,  77.144840,  79.337375,  82.910381,  84.735493,  87.425275,
    88.809111,  92.491899,  94.651344,  95.870634,  98.831194,  101.317851,
    103.736391, 105.446623, 107.168611, 110.434521, 111.472926, 113.747199,
    114.319730, 116.324265, 118.561636, 119.644908, 121.445660, 122.793985,
    124.578429, 125.744776, 127.556582, 129.653184, 130.727178, 132.610158,
    134.025335, 135.724500};

static const int PRIMES[] = {2,  3,  5,  7,  11, 13, 17, 19, 23, 29, 31, 37, 41,
                             43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
#define NPRIMES 25

static int deboor(const double *lam, const double *mu, int N, double *a,
                  double *b) {
  for (int k = 0; k < N - 1; k++)
    if (!(lam[k] < mu[k] && mu[k] < lam[k + 1]))
      return -1;
  double w[50], ws = 0.0;
  for (int k = 0; k < N; k++) {
    double num = 1.0;
    for (int j = 0; j < N - 1; j++)
      num *= lam[k] - mu[j];
    double den = 1.0;
    for (int j = 0; j < N; j++)
      if (j != k)
        den *= lam[k] - lam[j];
    w[k] = num / den;
    if (w[k] < 0.0)
      return -2;
    ws += w[k];
  }
  for (int k = 0; k < N; k++)
    w[k] /= ws;
  a[0] = 0.0;
  for (int i = 0; i < N; i++)
    a[0] += w[i] * lam[i];
  double np1 = 0.0;
  for (int i = 0; i < N; i++) {
    double v = lam[i] - a[0];
    np1 += w[i] * v * v;
  }
  b[0] = sqrt(np1);
  double npk = np1;
  for (int k = 1; k < N; k++) {
    double num = 0.0;
    for (int i = 0; i < N; i++) {
      double pp = 0.0, pc = 1.0;
      for (int j = 0; j < k; j++) {
        double b2 = (j > 0) ? b[j - 1] * b[j - 1] : 0.0;
        double pn = (lam[i] - a[j]) * pc - b2 * pp;
        pp = pc;
        pc = pn;
      }
      num += w[i] * lam[i] * pc * pc;
    }
    a[k] = num / npk;
    if (k < N - 1) {
      double npx = 0.0;
      for (int i = 0; i < N; i++) {
        double pp = 0.0, pc = 1.0;
        for (int j = 0; j <= k; j++) {
          double b2 = (j > 0) ? b[j - 1] * b[j - 1] : 0.0;
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

static int polyfit(const double *x, const double *y, int n, int deg,
                   double *c) {
  double ATA[36] = {0}, ATy[6] = {0};
  for (int i = 0; i < n; i++) {
    double xp[6];
    xp[0] = 1.0;
    for (int j = 1; j <= deg; j++)
      xp[j] = xp[j - 1] * x[i];
    for (int r = 0; r <= deg; r++) {
      ATy[r] += xp[r] * y[i];
      for (int s = 0; s <= deg; s++)
        ATA[r * (deg + 1) + s] += xp[r] * xp[s];
    }
  }
  for (int col = 0; col <= deg; col++) {
    int piv = col;
    double pmax = fabs(ATA[col * (deg + 1) + col]);
    for (int row = col + 1; row <= deg; row++) {
      double v = fabs(ATA[row * (deg + 1) + col]);
      if (v > pmax) {
        pmax = v;
        piv = row;
      }
    }
    if (pmax < 1e-30)
      return -1;
    if (piv != col) {
      for (int j = 0; j <= deg; j++) {
        double t = ATA[col * (deg + 1) + j];
        ATA[col * (deg + 1) + j] = ATA[piv * (deg + 1) + j];
        ATA[piv * (deg + 1) + j] = t;
      }
      {
        double t = ATy[col];
        ATy[col] = ATy[piv];
        ATy[piv] = t;
      }
    }
    double piv_val = ATA[col * (deg + 1) + col];
    for (int row = col + 1; row <= deg; row++) {
      double f = ATA[row * (deg + 1) + col] / piv_val;
      for (int j = col; j <= deg; j++)
        ATA[row * (deg + 1) + j] -= f * ATA[col * (deg + 1) + j];
      ATy[row] -= f * ATy[col];
    }
  }
  for (int i = deg; i >= 0; i--) {
    double s = ATy[i];
    for (int j = i + 1; j <= deg; j++)
      s -= ATA[i * (deg + 1) + j] * c[j];
    c[i] = s / ATA[i * (deg + 1) + i];
  }
  return 0;
}

static double polyval(const double *c, int deg, double x) {
  double y = c[deg];
  for (int i = deg - 1; i >= 0; i--)
    y = y * x + c[i];
  return y;
}

static void section(const char *s) { printf("\n=== %s ===\n\n", s); }

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

  printf(
      "==================================================================\n");
  printf("  Detrended Prime-Frequency Extraction (Euler form)\n");
  printf("  sin(k ln p) = [p^{ik} - p^{-ik}]/2i\n");
  printf("  DFT: F(omega) = (1/N) sum x_k * e^{-i*omega*k}\n");
  printf(
      "==================================================================\n\n");

  int N_test[] = {10, 15, 20, 25, 30, 35, 40, 45, 50};
  int n_test = 9;

  section("Method 1: Polynomial Detrend of Raw dBG Entries");

  for (int di = 3; di <= 5; di++) {
    printf("  --- Polynomial degree = %d ---\n\n", di);

    printf("  %6s  %12s  %12s  %12s  %12s\n", "N", "R2(a,sin)", "R2(b,cos)",
           "R2(comb)", "rms_resid");
    printf("  %6s  %12s  %12s  %12s  %12s\n", "---", "---", "---", "---",
           "---");

    for (int ti = 0; ti < n_test; ti++) {
      int N = N_test[ti];
      double lam[50], mu[50];
      for (int k = 0; k < N; k++)
        lam[k] = ZZ[k];
      for (int k = 0; k < N - 1; k++)
        mu[k] = 0.5 * (ZZ[k] + ZZ[k + 1]);
      double a[50], b[49];
      if (deboor(lam, mu, N, a, b) != 0)
        continue;

      double x[50];
      for (int k = 0; k < N; k++)
        x[k] = (double)k / (double)(N - 1);

      double ca[6] = {0}, cb[6] = {0};
      polyfit(x, a, N, di, ca);
      polyfit(x, b, N - 1, di, cb);

      double ra[50], rb[50];
      for (int k = 0; k < N; k++)
        ra[k] = a[k] - polyval(ca, di, x[k]);
      for (int k = 0; k < N - 1; k++)
        rb[k] = b[k] - polyval(cb, di, x[k]);

      int np = NPRIMES;
      if (np > N / 2)
        np = N / 2;

      double ssa = 0.0, sta = 0.0, ssb = 0.0, stb = 0.0;
      double ma = 0.0, mb = 0.0;
      for (int k = 0; k < N - 1; k++) {
        ma += ra[k];
        mb += rb[k];
      }
      ma /= (double)(N - 1);
      mb /= (double)(N - 1);

      for (int k = 0; k < N - 1; k++) {
        double pred_a = 0.0, pred_b = 0.0;
        for (int pi = 0; pi < np; pi++) {
          double p = (double)PRIMES[pi];
          double alpha = -log(p) / (2.0 * M_PI * sqrt(p));
          pred_a += alpha * euler_sin(p, (double)k);
          pred_b += alpha * euler_cos(p, (double)k);
        }
        ssa += (ra[k] - pred_a) * (ra[k] - pred_a);
        sta += (ra[k] - ma) * (ra[k] - ma);
        ssb += (rb[k] - pred_b) * (rb[k] - pred_b);
        stb += (rb[k] - mb) * (rb[k] - mb);
      }

      double r2a = (sta > 1e-30) ? 1.0 - ssa / sta : 0.0;
      double r2b = (stb > 1e-30) ? 1.0 - ssb / stb : 0.0;
      double r2c = (sta + stb > 1e-30) ? 1.0 - (ssa + ssb) / (sta + stb) : 0.0;

      double rms = sqrt((ssa + ssb) / (2.0 * (double)(N - 1)));

      printf("  %6d  %12.6f  %12.6f  %12.6f  %12.6f\n", N, r2a, r2b, r2c, rms);
    }
    printf("\n");
  }

  section("Method 2: First Differences (da_k, db_k)");

  printf("If a_k = smooth(k) + prime oscillation, then delta = high-freq "
         "part.\n\n");

  printf("  %6s  %12s  %12s  %12s  %12s\n", "N", "R2(da,sin)", "R2(db,cos)",
         "R2(comb)", "rms_resid");
  printf("  %6s  %12s  %12s  %12s  %12s\n", "---", "---", "---", "---", "---");

  for (int ti = 0; ti < n_test; ti++) {
    int N = N_test[ti];
    double lam[50], mu[50];
    for (int k = 0; k < N; k++)
      lam[k] = ZZ[k];
    for (int k = 0; k < N - 1; k++)
      mu[k] = 0.5 * (ZZ[k] + ZZ[k + 1]);
    double a[50], b[49];
    if (deboor(lam, mu, N, a, b) != 0)
      continue;

    int Nd = N - 2;
    double dda[50], ddb[50];
    for (int k = 0; k < Nd; k++) {
      dda[k] = a[k + 1] - a[k];
      ddb[k] = b[k + 1] - b[k];
    }

    int np = NPRIMES;
    if (np > Nd / 2)
      np = Nd / 2;

    double ssa = 0.0, sta = 0.0, ssb = 0.0, stb = 0.0;
    double ma = 0.0, mb = 0.0;
    for (int k = 0; k < Nd; k++) {
      ma += dda[k];
      mb += ddb[k];
    }
    ma /= (double)Nd;
    mb /= (double)Nd;

    for (int k = 0; k < Nd; k++) {
      double pred_a = 0.0, pred_b = 0.0;
      for (int pi = 0; pi < np; pi++) {
        double p = (double)PRIMES[pi];
        double alpha = -log(p) / (2.0 * M_PI * sqrt(p));
        double dw_sin = 2.0 * euler_cos(p, (double)k + 0.5) * euler_sin(p, 0.5);
        double dw_cos =
            -2.0 * euler_sin(p, (double)k + 0.5) * euler_sin(p, 0.5);
        pred_a += alpha * dw_sin;
        pred_b += alpha * dw_cos;
      }
      ssa += (dda[k] - pred_a) * (dda[k] - pred_a);
      sta += (dda[k] - ma) * (dda[k] - ma);
      ssb += (ddb[k] - pred_b) * (ddb[k] - pred_b);
      stb += (ddb[k] - mb) * (ddb[k] - mb);
    }

    double r2a = (sta > 1e-30) ? 1.0 - ssa / sta : 0.0;
    double r2b = (stb > 1e-30) ? 1.0 - ssb / stb : 0.0;
    double r2c = (sta + stb > 1e-30) ? 1.0 - (ssa + ssb) / (sta + stb) : 0.0;
    double rms = sqrt((ssa + ssb) / (2.0 * (double)Nd));

    printf("  %6d  %12.6f  %12.6f  %12.6f  %12.6f\n", N, r2a, r2b, r2c, rms);
  }

  section("Method 3: Per-Prime Coefficient Extraction (Poly Detrend deg=4)");

  int best_deg = 4;

  for (int ti = 4; ti < n_test; ti++) {
    int N = N_test[ti];
    double lam[50], mu[50];
    for (int k = 0; k < N; k++)
      lam[k] = ZZ[k];
    for (int k = 0; k < N - 1; k++)
      mu[k] = 0.5 * (ZZ[k] + ZZ[k + 1]);
    double a[50], b[49];
    if (deboor(lam, mu, N, a, b) != 0)
      continue;

    double x[50];
    for (int k = 0; k < N; k++)
      x[k] = (double)k / (double)(N - 1);

    double ca[6] = {0}, cb[6] = {0};
    polyfit(x, a, N, best_deg, ca);
    polyfit(x, b, N - 1, best_deg, cb);

    double ra[50], rb[50];
    for (int k = 0; k < N; k++)
      ra[k] = a[k] - polyval(ca, best_deg, x[k]);
    for (int k = 0; k < N - 1; k++)
      rb[k] = b[k] - polyval(cb, best_deg, x[k]);

    printf("  N = %d (poly detrend deg=%d):\n", N, best_deg);
    printf("  %4s  %12s  %12s  %12s  %12s  %8s  %12s  %12s\n", "p",
           "alpha_pred", "DFT_a.Re", "DFT_a.Im", "c_p(sin)", "ratio", "|DFT|",
           "conj_err");
    printf("  %4s  %12s  %12s  %12s  %12s  %8s  %12s  %12s\n", "---", "---",
           "---", "---", "---", "---", "---", "---");

    int np = NPRIMES;
    if (np > N / 2)
      np = N / 2;

    for (int pi = 0; pi < np; pi++) {
      double p = (double)PRIMES[pi];
      double omega = log(p);
      double alpha_pred = -log(p) / (2.0 * M_PI * sqrt(p));

      Cpx dft_a = dft_at_freq(ra, N - 1, omega);
      Cpx dft_b = dft_at_freq(rb, N - 1, omega);
      double cs_err = conj_symmetry_error(ra, N - 1, omega);

      double cp_sin = -2.0 * dft_a.im;
      double cp_cos_val = 2.0 * dft_b.re;
      double fmag = sqrt(cpx_abs2(dft_a) + cpx_abs2(dft_b));
      double ratio = (fabs(alpha_pred) > 1e-15) ? cp_sin / alpha_pred : 0.0;

      printf("  %4.0f  %12.8f  %12.8f  %12.8f  %12.8f  %8.4f  %12.8f  %12.2e  "
             "(cos=%.2e)\n",
             p, alpha_pred, dft_a.re, dft_a.im, cp_sin, ratio, fmag, cs_err,
             cp_cos_val);
    }
    printf("\n");
  }

  section("Method 4: Cross-N Residual Convergence at Fixed k");

  printf("For fixed k, detrend across N: residual(k,N) should converge.\n\n");

  double a_table[15][50], b_table[15][50];
  int valid[15] = {0};
  for (int ti = 0; ti < n_test; ti++) {
    int N = N_test[ti];
    double lam[50], mu[50];
    for (int k = 0; k < N; k++)
      lam[k] = ZZ[k];
    for (int k = 0; k < N - 1; k++)
      mu[k] = 0.5 * (ZZ[k] + ZZ[k + 1]);
    double a[50], b[49];
    if (deboor(lam, mu, N, a, b) != 0)
      continue;
    valid[ti] = 1;
    for (int k = 0; k < N; k++)
      a_table[ti][k] = a[k];
    for (int k = 0; k < N - 1; k++)
      b_table[ti][k] = b[k];
  }

  printf("  Detrended a_k residual (poly deg=4 in k) at fixed k:\n");
  printf("  %4s", "k\\N");
  for (int ti = 0; ti < n_test; ti++)
    if (valid[ti])
      printf("  %10d", N_test[ti]);
  printf("  %12s\n", "euler_pred");
  for (int k = 0; k < 6; k++) {
    printf("  %4d", k);
    for (int ti = 0; ti < n_test; ti++) {
      if (!valid[ti] || k >= N_test[ti]) {
        printf("  %10s", "---");
        continue;
      }
      int N = N_test[ti];
      double x[50];
      for (int j = 0; j < N; j++)
        x[j] = (double)j / (double)(N - 1);
      double ca[6] = {0};
      polyfit(x, a_table[ti], N, 4, ca);
      double res = a_table[ti][k] - polyval(ca, 4, x[k]);
      printf("  %10.6f", res);
    }
    double pred = 0.0;
    for (int pi = 0; pi < 12; pi++) {
      double p = (double)PRIMES[pi];
      double alpha = -log(p) / (2.0 * M_PI * sqrt(p));
      pred += alpha * euler_sin(p, (double)k);
    }
    printf("  %12.6f\n", pred);
  }

  printf("\n  Detrended b_k residual (poly deg=4 in k) at fixed k:\n");
  printf("  %4s", "k\\N");
  for (int ti = 0; ti < n_test; ti++)
    if (valid[ti])
      printf("  %10d", N_test[ti]);
  printf("  %12s\n", "euler_pred");
  for (int k = 0; k < 6; k++) {
    printf("  %4d", k);
    for (int ti = 0; ti < n_test; ti++) {
      if (!valid[ti] || k >= N_test[ti] - 1) {
        printf("  %10s", "---");
        continue;
      }
      int N = N_test[ti];
      double x[50];
      for (int j = 0; j < N - 1; j++)
        x[j] = (double)j / (double)(N - 2);
      double cb[6] = {0};
      polyfit(x, b_table[ti], N - 1, 4, cb);
      double res = b_table[ti][k] - polyval(cb, 4, x[k]);
      printf("  %10.6f", res);
    }
    double pred = 0.0;
    for (int pi = 0; pi < 12; pi++) {
      double p = (double)PRIMES[pi];
      double alpha = -log(p) / (2.0 * M_PI * sqrt(p));
      pred += alpha * euler_cos(p, (double)k);
    }
    printf("  %12.6f\n", pred);
  }

  section("Method 5: Least-Squares Coefficient Fit (N=50, poly detrend)");

  {
    int N = 50;
    double lam[50], mu[50];
    for (int k = 0; k < N; k++)
      lam[k] = ZZ[k];
    for (int k = 0; k < N - 1; k++)
      mu[k] = 0.5 * (ZZ[k] + ZZ[k + 1]);
    double a[50], b[49];
    if (deboor(lam, mu, N, a, b) != 0) {
      printf("  dBG failed\n");
      return 1;
    }

    double x[50];
    for (int k = 0; k < N; k++)
      x[k] = (double)k / (double)(N - 1);

    double ca[6], cb[6];
    polyfit(x, a, N, 4, ca);
    polyfit(x, b, N - 1, 4, cb);

    double ra[50], rb[50];
    for (int k = 0; k < N; k++)
      ra[k] = a[k] - polyval(ca, 4, x[k]);
    for (int k = 0; k < N - 1; k++)
      rb[k] = b[k] - polyval(cb, 4, x[k]);

    int np = 15;

    double ATA[225] = {0};
    double ATy[15] = {0};

    for (int pi = 0; pi < np; pi++) {
      double p = (double)PRIMES[pi];
      for (int pj = 0; pj < np; pj++) {
        double pj_val = (double)PRIMES[pj];
        double dot = 0.0;
        for (int k = 0; k < N - 1; k++) {
          double kd = (double)k;
          dot += euler_sin(p, kd) * euler_sin(pj_val, kd) +
                 euler_cos(p, kd) * euler_cos(pj_val, kd);
        }
        ATA[pi * np + pj] = dot;
      }
      double omega = log(p);
      Cpx dft_a = dft_at_freq(ra, N - 1, omega);
      Cpx dft_b = dft_at_freq(rb, N - 1, omega);
      double dot_y = -(double)(N - 1) * dft_a.im + (double)(N - 1) * dft_b.re;
      ATy[pi] = dot_y;
    }

    for (int col = 0; col < np; col++) {
      int piv = col;
      double pmax = fabs(ATA[col * np + col]);
      for (int row = col + 1; row < np; row++) {
        double v = fabs(ATA[row * np + col]);
        if (v > pmax) {
          pmax = v;
          piv = row;
        }
      }
      if (pmax < 1e-30)
        continue;
      if (piv != col) {
        for (int j = 0; j < np; j++) {
          double t = ATA[col * np + j];
          ATA[col * np + j] = ATA[piv * np + j];
          ATA[piv * np + j] = t;
        }
        {
          double t = ATy[col];
          ATy[col] = ATy[piv];
          ATy[piv] = t;
        }
      }
      double pv = ATA[col * np + col];
      for (int row = col + 1; row < np; row++) {
        double f = ATA[row * np + col] / pv;
        for (int j = col; j < np; j++)
          ATA[row * np + j] -= f * ATA[col * np + j];
        ATy[row] -= f * ATy[col];
      }
    }
    double c[15] = {0};
    for (int i = np - 1; i >= 0; i--) {
      double s = ATy[i];
      for (int j = i + 1; j < np; j++)
        s -= ATA[i * np + j] * c[j];
      if (fabs(ATA[i * np + i]) > 1e-30)
        c[i] = s / ATA[i * np + i];
    }

    printf("  Least-squares fit: find c_p minimizing |residual - sum "
           "c_p*f_p|^2\n");
    printf("  where f_p(k) = [euler_sin(p,k), euler_cos(p,k)]\n\n");

    printf("  %4s  %12s  %12s  %12s  %8s\n", "p", "alpha_pred", "c_p (LSQ)",
           "ratio", "sign");
    printf("  %4s  %12s  %12s  %12s  %8s\n", "---", "---", "---", "---", "---");

    double sum_sq_pred = 0.0, sum_sq_lsq = 0.0;
    for (int pi = 0; pi < np; pi++) {
      double p = (double)PRIMES[pi];
      double alpha_pred = -log(p) / (2.0 * M_PI * sqrt(p));
      double ratio = (fabs(alpha_pred) > 1e-15) ? c[pi] / alpha_pred : 0.0;
      int sign_ok = (c[pi] * alpha_pred > 0.0) ? 1 : 0;
      printf("  %4.0f  %12.8f  %12.8f  %12.4f  %8s\n", p, alpha_pred, c[pi],
             ratio, sign_ok ? "SAME" : "FLIP");
      sum_sq_pred += alpha_pred * alpha_pred;
      sum_sq_lsq += c[pi] * c[pi];
    }

    double ss_res = 0.0, ss_tot = 0.0;
    double mean_r = 0.0;
    for (int k = 0; k < N - 1; k++)
      mean_r += ra[k] + rb[k];
    mean_r /= (2.0 * (double)(N - 1));

    for (int k = 0; k < N - 1; k++) {
      double pred_a = 0.0, pred_b = 0.0;
      for (int pi = 0; pi < np; pi++) {
        double p = (double)PRIMES[pi];
        pred_a += c[pi] * euler_sin(p, (double)k);
        pred_b += c[pi] * euler_cos(p, (double)k);
      }
      ss_res += (ra[k] - pred_a) * (ra[k] - pred_a) +
                (rb[k] - pred_b) * (rb[k] - pred_b);
      ss_tot += (ra[k] - mean_r) * (ra[k] - mean_r) +
                (rb[k] - mean_r) * (rb[k] - mean_r);
    }
    double r2_lsq = (ss_tot > 1e-30) ? 1.0 - ss_res / ss_tot : 0.0;

    printf("\n  LSQ R^2 = %.6f\n", r2_lsq);
    printf("  ||c_lsq|| / ||alpha_pred|| = %.4f\n",
           sqrt(sum_sq_lsq / (sum_sq_pred + 1e-30)));
  }

  printf("\nDetrended analysis complete.\n");
  return 0;
}
