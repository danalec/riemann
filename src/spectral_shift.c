/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Spectral shift DFT: prime detection
 * @paper   yamaguchi-rh-2026.tex, Section 7.3
 * @theorem Theorem II
 * @proof   Prime-frequency encoding
 * @step    2
 *
 * spectral_shift.c — Spectral shift function xi(E) for midpoint-reconstructed
 *                    Jacobi at N=25. Fourier analysis at prime frequencies.
 *
 * xi(E) = N_zeta(E) - N_free(E)
 * where N_zeta(E) = Sturm count of reconstructed Jacobi eigenvalues <= E
 *       N_free(E) = theta(E)/pi   (Riemann-von Mangoldt smooth term)
 *
 * The Fourier transform of xi(E) should show peaks at omega = log p.
 * Explicit formula prediction: |F(log p)| ~ Delta_E / (2*pi*sqrt(p))
 *
 * Uses Euler exponential form:
 *   F(ω) = Σ xi(E_k) · e^{-iωE_k} · dE  (via dft_at_freq)
 *   sin(E ln p) = [p^{iE} - p^{-iE}] / 2i
 *
 * Compile: gcc -Wall -Wextra -Wconversion -Wshadow -Werror -O3
 *          -fno-strict-aliasing -fno-peel-loops -fno-unswitch-loops
 *          -I. -o spectral_shift.exe src/spectral_shift.c -lm
 */

#include "random_matrix_utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef _WIN32
#include <windows.h>
#endif

static const double zeta_zeros[30] = {
    14.134725141734695, 21.022039638771556, 25.010857580145689,
    30.424876125859512, 32.935061587739192, 37.586178158825675,
    40.918719012147498, 43.327073280915002, 48.005150881167161,
    49.773832477672300, 52.970321477714464, 56.446247697063392,
    59.347044002602352, 60.831778524609810, 65.112544048081602,
    67.079810529494168, 69.546401711173985, 72.067157674481905,
    75.704690699083926, 77.144840068874799, 79.337375020249368,
    82.910380854086029, 84.735492980517051, 87.425274613125225,
    88.809111207634459, 92.491899270558491, 94.651344040519888,
    95.870634228245308, 98.831194218193687, 101.317851005731384};

static double theta_s(double t) {
  if (t <= 1.0)
    return -M_PI / 8.0;
  double x = t / (2.0 * M_PI), u = 1.0 / t, u2 = u * u, u4 = u2 * u2;
  return 0.5 * t * log(x) - 0.5 * t - M_PI / 8.0 + u / 48.0 +
         7.0 * u * u2 / 5760.0 + 31.0 * u4 / 80640.0 +
         127.0 * u4 * u2 / 430080.0 + 2555.0 * u4 * u4 / 27525120.0 +
         1414477.0 * u4 * u4 * u2 / 18681062400.0;
}

static double theta_p(double t) {
  if (t <= 2.0 * M_PI)
    return 1.0;
  double t2 = t * t;
  return 0.5 * log(t / (2.0 * M_PI)) - 1.0 / (24.0 * t2) +
         7.0 / (960.0 * t2 * t2) + 31.0 / (8064.0 * t2 * t2 * t2);
}

static double gram(int n) {
  double g = (n == 0) ? 17.845599540410860
                      : 2.0 * M_PI * (double)n / log((double)n + 1.0);
  for (int i = 0; i < 30; i++) {
    double f = theta_s(g) - M_PI * (double)n;
    double d = f / theta_p(g);
    if (fabs(d) < 1e-15 * (1.0 + fabs(g)))
      break;
    g -= d;
  }
  return g;
}

static int deboor(const double *lam, const double *mu, int N, double *a,
                  double *b) {
  for (int k = 0; k < N - 1; k++)
    if (!(lam[k] < mu[k] && mu[k] < lam[k + 1]))
      return -1;

  double w[30] = {0}, ws = 0;
  for (int k = 0; k < N; k++) {
    double nm = 1;
    for (int j = 0; j < N - 1; j++)
      nm *= lam[k] - mu[j];
    double dn = 1;
    for (int j = 0; j < N; j++)
      if (j != k)
        dn *= lam[k] - lam[j];
    w[k] = nm / dn;
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
        double b2 = (j > 0) ? b[j - 1] * b[j - 1] : 0;
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
          double b2 = (j > 0) ? b[j - 1] * b[j - 1] : 0;
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

static int sturm_count(const double *d, const double *e, int N, double x) {
  int c = 0;
  double pp = 0.0, pc = 1.0;
  for (int k = 0; k < N; k++) {
    double ek = (k > 0) ? e[k - 1] : 0.0;
    double pn = (d[k] - x) * pc - ek * ek * pp;
    if (fabs(pn) > 1e100)
      pn = (pn > 0) ? 1e100 : -1e100;
    if (pc * pn < 0.0)
      c++;
    pp = pc;
    pc = pn;
  }
  return c;
}

static void jacobi_diag(int n, const double *a, const double *b, double *ev) {
  double *M = calloc((size_t)(n * n), sizeof(double));
  for (int k = 0; k < n; k++) {
    M[k * n + k] = a[k];
    if (k < n - 1)
      M[k * n + k + 1] = M[(k + 1) * n + k] = b[k];
  }
  for (int sw = 0; sw < 80; sw++) {
    double moff = 0;
    for (int p = 0; p < n - 1; p++)
      for (int q = p + 1; q < n; q++)
        if (fabs(M[p * n + q]) > moff)
          moff = fabs(M[p * n + q]);
    if (moff < 1e-15)
      break;
    for (int p = 0; p < n - 1; p++)
      for (int q = p + 1; q < n; q++) {
        double apq = M[p * n + q];
        double tol = 1e-16 * (fabs(M[p * n + p]) + fabs(M[q * n + q]) + 1);
        if (fabs(apq) < tol)
          continue;
        double app = M[p * n + p], aqq = M[q * n + q];
        double tau = (aqq - app) / (2.0 * apq);
        double t = (tau >= 0) ? 1.0 / (tau + sqrt(1.0 + tau * tau))
                              : -1.0 / (-tau + sqrt(1.0 + tau * tau));
        double c = 1.0 / sqrt(1.0 + t * t), s = t * c;
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
  free(M);
}

static int is_prime(int n) {
  if (n < 2)
    return 0;
  if (n == 2)
    return 1;
  if (n % 2 == 0)
    return 0;
  for (int d = 3; d * d <= n; d += 2)
    if (n % d == 0)
      return 0;
  return 1;
}

/* Trapezoidal DFT: F(ω) = dE · Σ_k xi_k · e^{-iωE_k} */
static Cpx trap_dft(const double *xi_arr, const double *E_arr, int n_samples,
                    double dE, double omega) {
  Cpx sum = cpx_make(0, 0);
  for (int i = 0; i < n_samples; i++) {
    Cpx twiddle = cpx_exp_i(-omega * E_arr[i]);
    sum = cpx_add(sum, cpx_scale(xi_arr[i], twiddle));
  }
  sum.re *= dE;
  sum.im *= dE;
  return sum;
}

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

  int N = 25;
  double lam[30], mu[29], a[30], b[29];

  for (int k = 0; k < N; k++)
    lam[k] = zeta_zeros[k];
  for (int k = 0; k < N - 1; k++)
    mu[k] = 0.5 * (zeta_zeros[k] + zeta_zeros[k + 1]);

  int ret = deboor(lam, mu, N, a, b);
  if (ret != 0) {
    fprintf(stderr, "de Boor-Golub failed (%d)\n", ret);
    return 1;
  }

  printf("=== Spectral Shift Function xi(E) = N_zeta(E) - theta(E)/pi ===\n");
  printf("=== Midpoint-reconstructed Jacobi, N = %d ===\n\n", N);
  printf("=== Using Euler exponential form: F(ω) = Σ xi·e^{-iωE}·dE ===\n\n");

  printf("De Boor-Golub reconstruction (a_k, b_k):\n");
  printf("  k   a_k          b_k\n");
  printf("  --- ----------  ----------\n");
  for (int k = 0; k < N; k++)
    printf("  %2d  %10.6f  %10.6f\n", k, a[k], (k < N - 1) ? b[k] : 0.0);

  double ev[30];
  jacobi_diag(N, a, b, ev);

  double max_err = 0;
  printf("\nEigenvalue verification:\n");
  printf("  k   gamma_k       lambda_k       error\n");
  printf("  --- ----------   ----------   ----------\n");
  for (int k = 0; k < N; k++) {
    double err = fabs(ev[k] - lam[k]);
    if (err > max_err)
      max_err = err;
    printf("  %2d  %10.6f   %10.6f   %.1e\n", k, lam[k], ev[k], err);
  }
  printf("  Max forward error: %.1e\n", max_err);

  /* ---- Compute xi(E) at 1000 sample points ---- */
  int n_samples = 1000;
  double E_min = lam[0] - 2.0;
  double E_max = lam[N - 1] + 5.0;
  double dE = (E_max - E_min) / (double)(n_samples - 1);

  double *xi_arr = malloc((size_t)n_samples * sizeof(double));
  double *E_arr = malloc((size_t)n_samples * sizeof(double));

  for (int i = 0; i < n_samples; i++) {
    E_arr[i] = E_min + (double)i * dE;
    int n_zeta = sturm_count(a, b, N, E_arr[i]);
    double n_free = theta_s(E_arr[i]) / M_PI;
    xi_arr[i] = (double)n_zeta - n_free;
  }

  printf("\nxi(E) sample values:\n");
  printf("  %8s  %10s  %10s  %10s\n", "E", "N_zeta", "theta/pi", "xi(E)");
  printf("  %8s  %10s  %10s  %10s\n", "--------", "--------", "--------",
         "--------");
  for (int i = 0; i < n_samples; i += 50) {
    int nz = sturm_count(a, b, N, E_arr[i]);
    double nf = theta_s(E_arr[i]) / M_PI;
    printf("  %8.3f  %10d  %10.4f  %10.4f\n", E_arr[i], nz, nf, xi_arr[i]);
  }

  printf("\nGram points g_n (first %d):\n", N);
  for (int n = 0; n < N; n++) {
    double g = gram(n);
    printf("  n=%2d  g_n=%10.6f  theta(g_n)/pi=%.1f\n", n, g,
           theta_s(g) / M_PI);
  }

  /*
   *  FOURIER TRANSFORM of xi(E) via Euler exponential form
   *  F(ω) = Σ_k xi(E_k) · e^{-iωE_k} · dE
   *
   *  For real xi(E), conjugate symmetry: F(-ω) = conj(F(+ω))
   */
  printf("\n==================================================================="
         "===\n");
  printf("  FOURIER SPECTRUM OF xi(E) AT PRIME FREQUENCIES (Euler form)\n");
  printf("====================================================================="
         "=\n\n");
  printf("  F(ω) = Σ xi(E_k) · e^{-iωE_k} · dE  (trapezoidal DFT)\n");
  printf("  Prediction: xi(E) ~ -(1/pi)*sum_p sin(E log p)/sqrt(p)\n");
  printf("  sin(E ln p) = [p^{iE} - p^{-iE}] / 2i\n");
  printf("  So |F(log p)| ~ Delta_E / (2*pi*sqrt(p)) = %.4f / sqrt(p)\n\n",
         (E_max - E_min) / (2.0 * M_PI));

  printf("  %5s  %10s  %12s  %12s  %12s  %12s  %12s  %10s\n", "p", "omega",
         "Re(F)", "Im(F)", "|F(omega)|", "predicted", "ratio", "peak?");
  printf("  %5s  %10s  %12s  %12s  %12s  %12s  %12s  %10s\n", "---", "---",
         "---", "---", "---", "---", "---", "---");

  double max_amp = 0;
  int max_p = 0;
  double delta_E = E_max - E_min;

  for (int p = 2; p <= 100; p++) {
    if (!is_prime(p))
      continue;

    double omega = log((double)p);
    Cpx F = trap_dft(xi_arr, E_arr, n_samples, dE, omega);
    double amp = cpx_abs(F);
    double pred = delta_E / (2.0 * M_PI * sqrt((double)p));
    double ratio = (pred > 1e-16) ? amp / pred : 0;

    if (amp > max_amp) {
      max_amp = amp;
      max_p = p;
    }

    printf("  %5d  %10.6f  %12.6f  %12.6f  %12.6f  %12.6f  %12.4f  %10s\n", p,
           omega, F.re, F.im, amp, pred, ratio, (amp > 0.05) ? "** PEAK" : "");
  }
  printf("\n  Max amplitude at p=%d: |F(log %d)| = %.6f\n", max_p, max_p,
         max_amp);

  /* Conjugate symmetry verification */
  printf("\n-------------------------------------------------------------------"
         "---\n");
  printf("  Conjugate symmetry check: F(-ω) should = conj(F(+ω))\n");
  printf("---------------------------------------------------------------------"
         "-\n\n");
  printf("  %5s  %12s  %12s  %12s  %12s\n", "p", "|F(+ω)|", "|F(-ω)|",
         "conj_err", "symmetric?");
  printf("  %5s  %12s  %12s  %12s  %12s\n", "---", "---", "---", "---", "---");
  {
    int check_primes[] = {2,  3,  5,  7,  11, 13, 17, 19, 23, 29, 31, 37, 41,
                          43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
    for (int i = 0; i < 25; i++) {
      double p = (double)check_primes[i];
      double omega = log(p);
      Cpx fp = trap_dft(xi_arr, E_arr, n_samples, dE, omega);
      Cpx fm = trap_dft(xi_arr, E_arr, n_samples, dE, -omega);
      Cpx expected_conj = cpx_conj(fp);
      double cerr = cpx_abs(cpx_sub(fm, expected_conj));
      printf("  %5d  %12.6f  %12.6f  %12.2e  %12s\n", check_primes[i],
             cpx_abs(fp), cpx_abs(fm), cerr, cerr < 0.01 ? "YES" : "NO");
    }
  }

  /* Comparison: prime vs non-prime */
  printf("\n-------------------------------------------------------------------"
         "---\n");
  printf("  Comparison: prime log-frequencies vs non-prime\n");
  printf("---------------------------------------------------------------------"
         "-\n\n");
  printf("  %12s  %12s  %12s  %12s  %8s\n", "omega", "Re(F)", "Im(F)",
         "|F(omega)|", "log(p)?");
  printf("  %12s  %12s  %12s  %12s  %8s\n", "---", "---", "---", "---", "---");

  double test_om[] = {log(2.0),  log(3.0),  log(5.0), log(7.0), log(11.0),
                      log(13.0), log(4.0),  log(6.0), log(8.0), log(9.0),
                      log(10.0), log(12.0), 1.0,      1.5,      2.0,
                      2.5,       3.0};
  int n_om = 17;

  for (int i = 0; i < n_om; i++) {
    double om = test_om[i];
    Cpx F = trap_dft(xi_arr, E_arr, n_samples, dE, om);

    int pv = (int)(exp(om) + 0.5);

    printf("  %12.6f  %12.6f  %12.6f  %12.6f  %8s\n", om, F.re, F.im,
           cpx_abs(F), is_prime(pv) ? "YES" : "no");
  }

  /* Summary: Fourier amplitudes at key primes */
  printf("\n==================================================================="
         "===\n");
  printf("  FOURIER AMPLITUDES |F(log p)| FOR SMALL PRIMES\n");
  printf("====================================================================="
         "=\n\n");

  int key_primes[] = {2,  3,  5,  7,  11, 13, 17, 19, 23, 29, 31, 37, 41,
                      43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
  int n_kp = 25;

  printf("  %3s  %12s  %12s  %12s  %10s\n", "p", "Re(F)", "Im(F)", "|F(log p)|",
         "ratio");
  printf("  %3s  %12s  %12s  %12s  %10s\n", "---", "---", "---", "---", "---");
  for (int i = 0; i < n_kp; i++) {
    double p = (double)key_primes[i];
    double om = log(p);
    Cpx F = trap_dft(xi_arr, E_arr, n_samples, dE, om);
    double pred = delta_E / (2.0 * M_PI * sqrt(p));
    printf("  %3d  %12.8f  %12.8f  %12.8f  %10.4f\n", key_primes[i], F.re, F.im,
           cpx_abs(F), cpx_abs(F) / pred);
  }

  free(xi_arr);
  free(E_arr);

  return 0;
}
