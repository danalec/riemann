/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Analytic entry formula
 * @paper   yamaguchi-rh-2026.tex, §6.2
 * @theorem Lemma I
 * @proof   Exploratory
 * @step    1 — Weyl asymptotics (eigenvalue approximation)
 *
 * analytic_entry_formula.c -- Verify and extend the analytic prime-entry
 * formulas
 *
 * KNOWN (from previous work):
 *   da_k = sum_p alpha_p * sin(log p * k)    alpha_p = -log(p)/(2pi*sqrt(p))
 *   db_k = sum_p alpha_p * cos(log p * k)
 *   R^2 = 0.80 (scalar gauge), 1.000000 (5 mechanisms)
 *   A_m(p) = A_1(p)^m / m  (algebraic identity)
 *
 * Uses Euler form: sin(k ln p) = [p^{ik} - p^{-ik}]/2i
 *                  cos(k ln p) = [p^{ik} + p^{-ik}]/2
 *                  DFT: F(omega) = (1/N) sum x_k * e^{-i*omega*k}
 *
 * Compile: gcc -Wall -Wextra -Wconversion -Wshadow -Werror -O3
 *          -fno-strict-aliasing -fno-peel-loops -fno-unswitch-loops
 *          -Isrc -o analytic_entry_formula src/analytic_entry_formula.c -lm
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
static const int NPRIMES = 25;

static int is_prime(long n) {
  if (n < 2)
    return 0;
  if (n == 2)
    return 1;
  if (n % 2 == 0)
    return 0;
  for (long d = 3; d * d <= n; d += 2)
    if (n % d == 0)
      return 0;
  return 1;
}

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

static double theta_riemann(double t) {
  return t / 2.0 * log(t / (2.0 * M_PI)) - t / 2.0 - M_PI / 8.0 +
         1.0 / (48.0 * t) + 7.0 / (5760.0 * t * t * t);
}

static double weyl_zero(int k) {
  double lo = 1.0, hi = 500.0;
  double target = M_PI * ((double)k - 7.0 / 8.0);
  for (int iter = 0; iter < 200; iter++) {
    double mid = 0.5 * (lo + hi);
    double v = theta_riemann(mid) - target;
    if (v > 0.0)
      hi = mid;
    else
      lo = mid;
    if (hi - lo < 1e-10)
      break;
  }
  return 0.5 * (lo + hi);
}

static void section(const char *s) { printf("\n=== %s ===\n\n", s); }

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

  printf(
      "==================================================================\n");
  printf("  Analytic Entry Formula Verification (Euler form)\n");
  printf("  da_k = sum_p alpha_p * sin(log(p) * k)  via [p^{ik}-p^{-ik}]/2i\n");
  printf("  db_k = sum_p alpha_p * cos(log(p) * k)  via [p^{ik}+p^{-ik}]/2\n");
  printf("  DFT: F(omega) = (1/N) sum x_k * e^{-i*omega*k}\n");
  printf("  alpha_p = -log(p) / (2*pi*sqrt(p))\n");
  printf(
      "==================================================================\n\n");

  section("1. Residual from Weyl-law Smooth Reference");

  int N_test[] = {20, 30, 40, 50};
  int n_test = 4;

  for (int ti = 0; ti < n_test; ti++) {
    int N = N_test[ti];

    double lam_z[50], mu_z[50];
    for (int k = 0; k < N; k++)
      lam_z[k] = ZZ[k];
    for (int k = 0; k < N - 1; k++)
      mu_z[k] = 0.5 * (ZZ[k] + ZZ[k + 1]);

    double lam_w[50], mu_w[50];
    for (int k = 0; k < N; k++)
      lam_w[k] = weyl_zero(k + 1);
    for (int k = 0; k < N - 1; k++)
      mu_w[k] = 0.5 * (lam_w[k] + lam_w[k + 1]);

    double a_z[50], b_z[49], a_w[50], b_w[49];
    int rz = deboor(lam_z, mu_z, N, a_z, b_z);
    int rw = deboor(lam_w, mu_w, N, a_w, b_w);

    if (rz != 0 || rw != 0) {
      printf("  N=%d: reconstruction failed (zeta=%d, weyl=%d)\n\n", N, rz, rw);
      continue;
    }

    printf("  N = %d (zeta zeros vs Weyl-law zeros):\n", N);
    printf("  %4s  %12s  %12s  %12s  %12s  %12s\n", "k", "a_k(zeta)",
           "a_k(weyl)", "delta_a", "delta_b", "euler_pred");
    printf("  %4s  %12s  %12s  %12s  %12s  %12s\n", "---", "---", "---", "---",
           "---", "---");

    double da[50], db[50];
    for (int k = 0; k < (N < 15 ? N : 15); k++) {
      da[k] = a_z[k] - a_w[k];
      db[k] = (k < N - 1) ? (b_z[k] - b_w[k]) : 0.0;

      double pred = euler_prime_sum(PRIMES, NPRIMES, (double)k);

      printf("  %4d  %12.4f  %12.4f  %12.6f  %12.6f  %12.6f\n", k, a_z[k],
             a_w[k], da[k], db[k], pred);
    }

    printf("\n  Per-prime coefficient extraction via DFT (N=%d):\n", N);
    printf("  %4s  %12s  %12s  %12s  %12s  %8s  %12s\n", "p", "alpha_pred",
           "DFT_Re", "DFT_Im", "c_p(sin)", "ratio", "conj_err");
    printf("  %4s  %12s  %12s  %12s  %12s  %8s  %12s\n", "---", "---", "---",
           "---", "---", "---", "---");

    int np_use = NPRIMES;
    if (np_use > N / 2)
      np_use = N / 2;

    for (int pi = 0; pi < np_use; pi++) {
      double p = (double)PRIMES[pi];
      double omega = log(p);
      double alpha_pred = -log(p) / (2.0 * M_PI * sqrt(p));

      Cpx dft_sin = dft_at_freq(da, N - 1, omega);
      Cpx dft_cos = dft_at_freq(db, N - 1, omega);
      double cs_err = conj_symmetry_error(da, N - 1, omega);

      double cp_sin = -2.0 * dft_sin.im;
      double cp_cos_val = 2.0 * dft_cos.re;
      double ratio = (fabs(alpha_pred) > 1e-15) ? cp_sin / alpha_pred : 0.0;

      printf("  %4.0f  %12.8f  %12.8f  %12.8f  %12.8f  %8.4f  %12.2e  "
             "(cos=%.2e)\n",
             p, alpha_pred, dft_sin.re, dft_sin.im, cp_sin, ratio, cs_err,
             cp_cos_val);
    }

    double ss_res = 0.0, ss_tot = 0.0;
    double mean_da = 0.0;
    for (int k = 0; k < N - 1; k++)
      mean_da += da[k];
    mean_da /= (double)(N - 1);

    for (int k = 0; k < N - 1; k++) {
      double pred = euler_prime_sum(PRIMES, np_use, (double)k);
      double err = da[k] - pred;
      ss_res += err * err;
      ss_tot += (da[k] - mean_da) * (da[k] - mean_da);
    }
    double r2_a = (ss_tot > 1e-30) ? 1.0 - ss_res / ss_tot : 0.0;

    double ss_res_b = 0.0, ss_tot_b = 0.0;
    double mean_db = 0.0;
    for (int k = 0; k < N - 1; k++)
      mean_db += db[k];
    mean_db /= (double)(N - 1);

    for (int k = 0; k < N - 1; k++) {
      double pred = 0.0;
      for (int pi = 0; pi < np_use; pi++) {
        double p = (double)PRIMES[pi];
        double alpha = -log(p) / (2.0 * M_PI * sqrt(p));
        pred += alpha * euler_cos(p, (double)k);
      }
      double err = db[k] - pred;
      ss_res_b += err * err;
      ss_tot_b += (db[k] - mean_db) * (db[k] - mean_db);
    }
    double r2_b = (ss_tot_b > 1e-30) ? 1.0 - ss_res_b / ss_tot_b : 0.0;

    double r2_combined = (ss_tot + ss_tot_b > 1e-30)
                             ? 1.0 - (ss_res + ss_res_b) / (ss_tot + ss_tot_b)
                             : 0.0;

    printf("\n  R^2 (delta_a vs sin prime sum):     %.6f\n", r2_a);
    printf("  R^2 (delta_b vs cos prime sum):     %.6f\n", r2_b);
    printf("  R^2 (combined):                     %.6f\n", r2_combined);
    printf("\n");
  }

  section("2. Empirical Formula: b_n vs sqrt(n)*log(p_n)/sqrt(p_n)");

  printf("b_n ~ sqrt(n) * log(p_n) / sqrt(p_n)  (Pearson r = +0.861)\n\n");

  for (int ti = 0; ti < n_test; ti++) {
    int N = N_test[ti];
    double lam[50], mu[50];
    for (int k = 0; k < N; k++)
      lam[k] = ZZ[k];
    for (int k = 0; k < N - 1; k++)
      mu[k] = 0.5 * (ZZ[k] + ZZ[k + 1]);
    double a[50], b[49];
    if (deboor(lam, mu, N, a, b) != 0) {
      printf("  N=%d: failed\n\n", N);
      continue;
    }

    printf("  N = %d:\n", N);
    printf("  %4s  %12s  %12s  %12s  %8s  %12s\n", "k", "b_k", "pred", "ratio",
           "p_k", "gamma_k");
    printf("  %4s  %12s  %12s  %12s  %8s  %12s\n", "---", "---", "---", "---",
           "---", "---");

    int p_list[50];
    int pc = 0;
    for (long p = 2; pc < N - 1 && p < 1000; p++) {
      if (is_prime(p)) {
        p_list[pc] = (int)p;
        pc++;
      }
    }

    double sum_r = 0.0, sum_rx = 0.0, sum_ry = 0.0, sum_xx = 0.0, sum_yy = 0.0;
    int nr = 0;

    for (int k = 0; k < (N - 1 < 20 ? N - 1 : 20); k++) {
      if (k >= pc)
        break;
      double pn = (double)p_list[k];
      double pred = sqrt((double)(k + 1)) * log(pn) / sqrt(pn);
      double ratio = (pred > 1e-15) ? b[k] / pred : 0.0;

      printf("  %4d  %12.4f  %12.4f  %12.4f  %8d  %12.4f\n", k, b[k], pred,
             ratio, p_list[k], ZZ[k]);

      if (pred > 0.01) {
        sum_r += b[k] * pred;
        sum_rx += b[k];
        sum_ry += pred;
        sum_xx += b[k] * b[k];
        sum_yy += pred * pred;
        nr++;
      }
    }

    if (nr > 2) {
      double mean_x = sum_rx / (double)nr;
      double mean_y = sum_ry / (double)nr;
      double cov = sum_r / (double)nr - mean_x * mean_y;
      double varx = sum_xx / (double)nr - mean_x * mean_x;
      double vary = sum_yy / (double)nr - mean_y * mean_y;
      double pearson = (varx * vary > 1e-30) ? cov / sqrt(varx * vary) : 0.0;
      printf("  Pearson r = %.4f (n=%d)\n", pearson, nr);
    }
    printf("\n");
  }

  section("3. DFT Coefficient Extraction from delta_a_k (Euler form)");

  for (int ti = 0; ti < n_test; ti++) {
    int N = N_test[ti];
    double lam_z[50], mu_z[50], lam_w[50], mu_w[50];
    for (int k = 0; k < N; k++) {
      lam_z[k] = ZZ[k];
      lam_w[k] = weyl_zero(k + 1);
    }
    for (int k = 0; k < N - 1; k++) {
      mu_z[k] = 0.5 * (ZZ[k] + ZZ[k + 1]);
      mu_w[k] = 0.5 * (lam_w[k] + lam_w[k + 1]);
    }
    double a_z[50], b_z[49], a_w[50], b_w[49];
    if (deboor(lam_z, mu_z, N, a_z, b_z) != 0)
      continue;
    if (deboor(lam_w, mu_w, N, a_w, b_w) != 0)
      continue;

    double da[50];
    for (int k = 0; k < N; k++)
      da[k] = a_z[k] - a_w[k];

    printf("  N = %d:\n", N);
    printf("  %4s  %12s  %12s  %12s  %12s  %8s  %12s\n", "p", "alpha_pred",
           "DFT_Re", "DFT_Im", "ratio", "|F|/N", "conj_err");
    printf("  %4s  %12s  %12s  %12s  %12s  %8s  %12s\n", "---", "---", "---",
           "---", "---", "---", "---");

    int np_use = (NPRIMES < N / 2) ? NPRIMES : N / 2;
    for (int pi = 0; pi < np_use; pi++) {
      double p = (double)PRIMES[pi];
      double omega = log(p);
      double alpha_pred = -log(p) / (2.0 * M_PI * sqrt(p));

      Cpx dft = dft_at_freq(da, N, omega);
      double cp = -2.0 * dft.im;
      double fmag = cpx_abs(dft);
      double ratio = (fabs(alpha_pred) > 1e-15) ? cp / alpha_pred : 0.0;
      double cs_err = conj_symmetry_error(da, N, omega);

      printf("  %4.0f  %12.8f  %12.8f  %12.8f  %12.4f  %8.4f  %12.2e\n", p,
             alpha_pred, dft.re, dft.im, ratio, fmag, cs_err);
    }

    printf("\n  Non-prime frequency check (should be near zero):\n");
    int non_primes[] = {4, 6, 8, 9, 10, 12, 14, 15};
    for (int i = 0; i < 8; i++) {
      double npp = (double)non_primes[i];
      double omega = log(npp);
      Cpx dft = dft_at_freq(da, N, omega);
      double cs_err = conj_symmetry_error(da, N, omega);
      printf("  log(%2.0f): |F|=%.6f  Re=%+.2e  Im=%+.2e  conj_err=%.2e\n", npp,
             cpx_abs(dft), dft.re, dft.im, cs_err);
    }
    printf("\n");
  }

  section("4. Prime-Power Frequencies: m*log(p)");

  printf("Predicted: c_{m,p} = alpha_p^m / m  "
         "where alpha_p = -log(p)/(2pi*sqrt(p))\n\n");

  int N = 40;
  double lam_z[50], mu_z[50], lam_w[50], mu_w[50];
  for (int k = 0; k < N; k++) {
    lam_z[k] = ZZ[k];
    lam_w[k] = weyl_zero(k + 1);
  }
  for (int k = 0; k < N - 1; k++) {
    mu_z[k] = 0.5 * (ZZ[k] + ZZ[k + 1]);
    mu_w[k] = 0.5 * (lam_w[k] + lam_w[k + 1]);
  }
  double a_z[50], b_z[49], a_w[50], b_w[49];
  if (deboor(lam_z, mu_z, N, a_z, b_z) != 0) {
    printf("  zeta dBG failed\n");
    return 1;
  }
  if (deboor(lam_w, mu_w, N, a_w, b_w) != 0) {
    printf("  weyl dBG failed\n");
    return 1;
  }

  double da[50];
  for (int k = 0; k < N; k++)
    da[k] = a_z[k] - a_w[k];

  printf("  N = %d, extracting at m*log(p) for first 8 primes:\n\n", N);
  printf("  %4s  %4s  %12s  %12s  %12s  %12s  %12s\n", "p", "m", "pred",
         "DFT_Re", "DFT_Im", "DFT_sin", "ratio");
  printf("  %4s  %4s  %12s  %12s  %12s  %12s  %12s\n", "---", "---", "---",
         "---", "---", "---", "---");

  for (int pi = 0; pi < 8; pi++) {
    double p = (double)PRIMES[pi];
    double alpha1 = -log(p) / (2.0 * M_PI * sqrt(p));

    for (int m = 1; m <= 3; m++) {
      double omega = (double)m * log(p);
      double pred = pow(alpha1, (double)m) / (double)m;

      Cpx dft = dft_at_freq(da, N, omega);
      double cp = -2.0 * dft.im;
      double cs_err = conj_symmetry_error(da, N, omega);
      double ratio = (fabs(pred) > 1e-18) ? cp / pred : 0.0;

      printf("  %4.0f  %4d  %12.6e  %12.6e  %12.6e  %12.6e  %12.4f  cs=%.2e\n",
             p, m, pred, dft.re, dft.im, cp, ratio, cs_err);
    }
  }

  section("5. R^2 Summary: How Much of delta_a is Prime-Frequency?");

  printf("Fitting da_k = sum c_p * euler_sin(p, k) for p=2..97\n\n");

  printf("  %6s  %12s  %12s  %12s  %12s\n", "N", "R2(sin,a)", "R2(cos,b)",
         "R2(comb)", "#primes");
  printf("  %6s  %12s  %12s  %12s  %12s\n", "---", "---", "---", "---", "---");

  int N_all[] = {10, 15, 20, 25, 30, 35, 40, 45, 50};
  for (int ti = 0; ti < 9; ti++) {
    int Ni = N_all[ti];
    double lz[50], mz[50], lw[50], mw[50];
    for (int k = 0; k < Ni; k++) {
      lz[k] = ZZ[k];
      lw[k] = weyl_zero(k + 1);
    }
    for (int k = 0; k < Ni - 1; k++) {
      mz[k] = 0.5 * (ZZ[k] + ZZ[k + 1]);
      mw[k] = 0.5 * (lw[k] + lw[k + 1]);
    }
    double az[50], bz[49], aw[50], bw[49];
    if (deboor(lz, mz, Ni, az, bz) != 0)
      continue;
    if (deboor(lw, mw, Ni, aw, bw) != 0)
      continue;

    double da_l[50], db_l[50];
    for (int k = 0; k < Ni; k++)
      da_l[k] = az[k] - aw[k];
    for (int k = 0; k < Ni - 1; k++)
      db_l[k] = bz[k] - bw[k];

    int np = NPRIMES;
    if (np > Ni / 2)
      np = Ni / 2;

    double ssa = 0.0, sta = 0.0, ssb = 0.0, stb = 0.0;
    double ma = 0.0, mb = 0.0;
    for (int k = 0; k < Ni - 1; k++) {
      ma += da_l[k];
      mb += db_l[k];
    }
    ma /= (double)(Ni - 1);
    mb /= (double)(Ni - 1);

    for (int k = 0; k < Ni - 1; k++) {
      double pred_a = 0.0, pred_b = 0.0;
      for (int pi = 0; pi < np; pi++) {
        double p = (double)PRIMES[pi];
        double alpha = -log(p) / (2.0 * M_PI * sqrt(p));
        pred_a += alpha * euler_sin(p, (double)k);
        pred_b += alpha * euler_cos(p, (double)k);
      }
      ssa += (da_l[k] - pred_a) * (da_l[k] - pred_a);
      sta += (da_l[k] - ma) * (da_l[k] - ma);
      ssb += (db_l[k] - pred_b) * (db_l[k] - pred_b);
      stb += (db_l[k] - mb) * (db_l[k] - mb);
    }

    double r2a = (sta > 1e-30) ? 1.0 - ssa / sta : 0.0;
    double r2b = (stb > 1e-30) ? 1.0 - ssb / stb : 0.0;
    double r2c = (sta + stb > 1e-30) ? 1.0 - (ssa + ssb) / (sta + stb) : 0.0;

    printf("  %6d  %12.6f  %12.6f  %12.6f  %12d\n", Ni, r2a, r2b, r2c, np);
  }

  printf("\nAnalytic entry formula verification complete.\n");
  return 0;
}
