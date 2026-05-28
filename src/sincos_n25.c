/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Sin/cos channel decomposition
 * @paper   yamaguchi-rh-2026.tex, Section 7.3
 * @theorem Theorem II
 * @proof   Fejer DFT
 * @step    2
 *
 * sincos_n25.c -- sin/cos first-order perturbation at N=25
 * Reconstructs Jacobi matrix via de Boor-Golub (midpoints),
 * computes eigenvectors via Jacobi with accumulation,
 * evaluates diagonal (sin) and off-diagonal (cos) channel strength.
 *
 * Uses Euler form: sin(k ln p) = [p^{ik} - p^{-ik}]/2i
 *                  cos(k ln p) = [p^{ik} + p^{-ik}]/2
 *
 * Compile: gcc -O3 -fno-strict-aliasing -fno-peel-loops -fno-unswitch-loops
 *               -Wall -Wextra -Wconversion -Wshadow -Werror -Isrc -lm
 *               -o sincos_n25 src/sincos_n25.c
 */

#include "random_matrix_utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef _WIN32
#include <windows.h>
#endif

#define N 25
#define PRIME_COUNT 12

static const double zeta[30] = {
    14.134725141734693790, 21.022039638771554993, 25.010857580145688763,
    30.424876125859513210, 32.935061587739189691, 37.586178158825671257,
    40.918719012147495187, 43.327073280914999519, 48.005150881167159728,
    49.773832477672302182, 52.970321477714460644, 56.446247697063394804,
    59.347044002602353080, 60.831778524609809844, 65.112544048081606661,
    67.079810529494173714, 69.546401711173979253, 72.067157674481907600,
    75.704690699083933168, 77.144840068874805373, 79.337375020249367923,
    82.910380854086030183, 84.735492980517050106, 87.425274613125229407,
    88.809111207634465424, 92.491899270558484295, 94.651344040519886967,
    95.870634228245309187, 98.831194218193692235, 101.317851005731391326};

static const int primes[PRIME_COUNT] = {2,  3,  5,  7,  11, 13,
                                        17, 19, 23, 29, 31, 37};

static double gram_theta(double t) {
  return (t / 2.0) * log(t / (2.0 * M_PI)) - t / 2.0 - M_PI / 8.0 +
         1.0 / (48.0 * t) + 7.0 / (5760.0 * t * t * t);
}

static double gram_point(int n) {
  double t0 = 2.0 * M_PI * exp(1.0);
  for (int it = 0; it < 30; it++) {
    double th = gram_theta(t0);
    double dth = 0.5 * log(t0 / (2.0 * M_PI)) - 1.0 / (48.0 * t0 * t0) -
                 21.0 / (5760.0 * t0 * t0 * t0 * t0);
    double f = th - (double)n * M_PI;
    double dt = f / dth;
    t0 -= dt;
    if (fabs(dt) < 1e-12)
      break;
  }
  return t0;
}

static int deboor(const double *lam, const double *mu, int n, double *a,
                  double *b) {
  for (int k = 0; k < n - 1; k++)
    if (!(lam[k] < mu[k] && mu[k] < lam[k + 1]))
      return -1;
  double w[30], ws = 0.0;
  for (int k = 0; k < n; k++) {
    double num = 1.0;
    for (int j = 0; j < n - 1; j++)
      num *= lam[k] - mu[j];
    double den = 1.0;
    for (int j = 0; j < n; j++)
      if (j != k)
        den *= lam[k] - lam[j];
    w[k] = num / den;
    if (w[k] < 0.0)
      return -2;
    ws += w[k];
  }
  for (int k = 0; k < n; k++)
    w[k] /= ws;
  a[0] = 0.0;
  for (int i = 0; i < n; i++)
    a[0] += w[i] * lam[i];
  double np1 = 0.0;
  for (int i = 0; i < n; i++) {
    double v = lam[i] - a[0];
    np1 += w[i] * v * v;
  }
  b[0] = sqrt(np1);
  double npk = np1;
  for (int k = 1; k < n; k++) {
    double num = 0.0;
    for (int i = 0; i < n; i++) {
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
    if (k < n - 1) {
      double npx = 0.0;
      for (int i = 0; i < n; i++) {
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

static void jacobi_eig(int n, const double *a, const double *b, double *evals,
                       double *U) {
  double *M = (double *)xcalloc((size_t)(n * n), sizeof(double));
  for (int i = 0; i < n; i++) {
    M[i * n + i] = a[i];
    if (i < n - 1) {
      M[i * n + (i + 1)] = b[i];
      M[(i + 1) * n + i] = b[i];
    }
  }
  for (int i = 0; i < n * n; i++)
    U[i] = 0.0;
  for (int i = 0; i < n; i++)
    U[i * n + i] = 1.0;

  for (int sw = 0; sw < 60; sw++) {
    double moff = 0.0;
    for (int p = 0; p < n - 1; p++)
      for (int q = p + 1; q < n; q++) {
        double v = fabs(M[p * n + q]);
        if (v > moff)
          moff = v;
      }
    if (moff < 1e-14)
      break;
    for (int p = 0; p < n - 1; p++) {
      for (int q = p + 1; q < n; q++) {
        double apq = M[p * n + q];
        if (fabs(apq) < 1e-16 * (fabs(M[p * n + p]) + fabs(M[q * n + q]) + 1.0))
          continue;
        double app = M[p * n + p], aqq = M[q * n + q];
        double tau = (aqq - app) / (2.0 * apq);
        double t = (tau >= 0.0) ? 1.0 / (tau + sqrt(1.0 + tau * tau))
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
        for (int i = 0; i < n; i++) {
          double uip = U[i * n + p], uiq = U[i * n + q];
          U[i * n + p] = uip * c - uiq * s;
          U[i * n + q] = uip * s + uiq * c;
        }
      }
    }
  }
  for (int i = 0; i < n; i++)
    evals[i] = M[i * n + i];
  for (int i = 1; i < n; i++) {
    double key = evals[i];
    int j = i - 1;
    double *col = (double *)xmalloc((size_t)n * sizeof(double));
    for (int r = 0; r < n; r++)
      col[r] = U[r * n + i];
    while (j >= 0 && evals[j] > key) {
      evals[j + 1] = evals[j];
      for (int r = 0; r < n; r++)
        U[r * n + (j + 1)] = U[r * n + j];
      j--;
    }
    evals[j + 1] = key;
    for (int r = 0; r < n; r++)
      U[r * n + (j + 1)] = col[r];
    free(col);
  }
  free(M);
}

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

  double lam[30], mu[29], a[30], b[29];
  for (int k = 0; k < N; k++)
    lam[k] = zeta[k];
  for (int k = 0; k < N - 1; k++)
    mu[k] = 0.5 * (zeta[k] + zeta[k + 1]);

  int ret = deboor(lam, mu, N, a, b);
  if (ret != 0) {
    fprintf(stderr, "deboor failed: %d\n", ret);
    return 1;
  }

  double evals[30], U[900];
  jacobi_eig(N, a, b, evals, U);

  double max_err = 0.0;
  for (int k = 0; k < N; k++) {
    double err = fabs(evals[k] - lam[k]);
    if (err > max_err)
      max_err = err;
  }
  printf("de Boor-Golub reconstruction error (max |ev - zeta|): %.2e\n\n",
         max_err);

  double alpha[PRIME_COUNT];
  for (int i = 0; i < PRIME_COUNT; i++)
    alpha[i] = -log((double)primes[i]) / (2.0 * M_PI * sqrt((double)primes[i]));

  double da[30], db[29];
  for (int j = 0; j < N; j++) {
    da[j] = 0.0;
    for (int i = 0; i < PRIME_COUNT; i++)
      da[j] += alpha[i] * euler_sin((double)primes[i], (double)j);
  }
  for (int j = 0; j < N - 1; j++) {
    db[j] = 0.0;
    for (int i = 0; i < PRIME_COUNT; i++)
      db[j] += alpha[i] * euler_cos((double)primes[i], (double)j);
  }

  double sin_sum[30], cos_sum[30];
  int s_count = 0, c_count = 0;

  printf(
      " k   eigenvalue      sin_sum       cos_sum        channel  |sin/cos|\n");
  printf(
      "---  -------------  ------------  ------------   -------  ---------\n");

  for (int k = 0; k < N; k++) {
    sin_sum[k] = 0.0;
    for (int j = 0; j < N; j++) {
      double u2 = U[j * N + k] * U[j * N + k];
      sin_sum[k] += u2 * da[j];
    }
    cos_sum[k] = 0.0;
    for (int j = 0; j < N - 1; j++) {
      cos_sum[k] += 2.0 * U[j * N + k] * U[(j + 1) * N + k] * db[j];
    }
    char channel = (fabs(sin_sum[k]) > fabs(cos_sum[k])) ? 's' : 'c';
    if (channel == 's')
      s_count++;
    else
      c_count++;
    double ratio = fabs(sin_sum[k]) / (fabs(cos_sum[k]) + 1e-30);
    printf("%2d   %12.6f   %+.7e   %+.7e     %c       %.3f\n", k, evals[k],
           sin_sum[k], cos_sum[k], channel, ratio);
  }

  double frac_s = (double)s_count / (double)N;
  double frac_c = (double)c_count / (double)N;

  printf("\n--- Summary ---\n");
  printf("N = %d\n", N);
  printf("sin-dominant eigenvalues: %d / %d = %.1f%%\n", s_count, N,
         100.0 * frac_s);
  printf("cos-dominant eigenvalues: %d / %d = %.1f%%\n", c_count, N,
         100.0 * frac_c);

  printf("\n--- Site-by-site perturbation (da_j, db_j) ---\n");
  printf(" j   site (midpt)    da_j (sin)    db_j (cos)\n");
  printf("---  --------------  ------------  ------------\n");
  for (int j = 0; j < N; j++) {
    printf("%2d   %12.6f   %+.7e   %+.7e\n", j + 1, mu[j], da[j],
           (j < N - 1) ? db[j] : 0.0);
  }

  printf("\n--- Fundamental trig identity: da_p[j]^2 + db_p[j]^2 = alpha_p^2 "
         "---\n");
  printf("  p    omega_p      alpha_p     alpha_p^2     max|err|     "
         "site(worst)\n");
  printf(
      "---  ---------   ----------   ----------   ----------   -----------\n");
  for (int pi = 0; pi < PRIME_COUNT; pi++) {
    double pval = (double)primes[pi];
    double ap = alpha[pi];
    double ap2 = ap * ap;
    double max_err_p = 0.0;
    int worst_j = 0;
    for (int j = 0; j < N; j++) {
      double da_p = ap * euler_sin(pval, (double)j);
      double db_p = (j < N - 1) ? ap * euler_cos(pval, (double)j) : 0.0;
      double sum2 = da_p * da_p + db_p * db_p;
      double err = fabs(sum2 - ap2);
      if (err > max_err_p) {
        max_err_p = err;
        worst_j = j;
      }
    }
    printf("%3d  %9.4f   %+10.6f   %10.6f   %10.2e   %3d\n", primes[pi],
           log(pval), ap, ap2, max_err_p, worst_j);
  }

  double sum_alpha2 = 0.0;
  for (int pi = 0; pi < PRIME_COUNT; pi++)
    sum_alpha2 += alpha[pi] * alpha[pi];

  printf("\n--- Per-site power da[j]^2+db[j]^2 vs sum(alpha_p^2) ---\n");
  printf("  j    da[j]^2+db[j]^2   sum_ap2     ratio    deviation\n");
  printf("---   ---------------  --------   ------   ---------\n");
  for (int j = 0; j < N; j++) {
    double pow_j = da[j] * da[j] + ((j < N - 1) ? db[j] * db[j] : 0.0);
    double dev = pow_j - sum_alpha2;
    printf("%3d   %14.6f  %9.6f   %6.3f   %+.6f\n", j, pow_j, sum_alpha2,
           pow_j / sum_alpha2, dev);
  }

  printf("\n--- Angle theta_j = atan2(da[j], db[j]) ---\n");
  printf("  j    theta_j      theta_j/pi   expected_2  expected_3  match?\n");
  printf("---   --------     ----------   ----------  ----------  ------\n");
  for (int j = 0; j < N && j < 15; j++) {
    double dbj = (j < N - 1) ? db[j] : 0.0;
    double theta_j = atan2(da[j], dbj);
    double omega2 = log(2.0), omega3 = log(3.0);
    double exp2 = fmod(omega2 * (double)j, 2.0 * M_PI);
    double exp3 = fmod(omega3 * (double)j, 2.0 * M_PI);
    if (exp2 < 0.0)
      exp2 += 2.0 * M_PI;
    if (exp3 < 0.0)
      exp3 += 2.0 * M_PI;
    double th_pos = theta_j;
    if (th_pos < 0.0)
      th_pos += 2.0 * M_PI;
    int m2 =
        (fabs(th_pos - exp2) < 0.3 || fabs(th_pos - exp2 - 2.0 * M_PI) < 0.3);
    int m3 =
        (fabs(th_pos - exp3) < 0.3 || fabs(th_pos - exp3 - 2.0 * M_PI) < 0.3);
    printf("%3d   %8.4f     %8.4f    %8.4f    %8.4f    %s%s\n", j, theta_j,
           theta_j / M_PI, exp2, exp3, m2 ? "p2" : "  ", m3 ? "p3" : "  ");
  }

  printf("\n--- Eigenvalue trig identity: sin_sum^2 + cos_sum^2 ---\n");
  printf("  k   sin^2+cos^2     sqrt       per-ev power\n");
  printf("---   -----------   --------   -------------\n");
  double mean_pow = 0.0;
  for (int k = 0; k < N; k++) {
    double pow_k = sin_sum[k] * sin_sum[k] + cos_sum[k] * cos_sum[k];
    mean_pow += pow_k;
    printf("%3d   %11.6f   %8.6f\n", k, pow_k, sqrt(pow_k));
  }
  mean_pow /= (double)N;
  printf("Mean power per eigenvalue: %.6f\n", mean_pow);

  printf("\n--- DFT at prime frequencies from perturbation vector da ---\n");
  printf("  p     omega       DFT_Re       DFT_Im      |DFT|      conj_err\n");
  printf("---  ---------   ----------   ----------   ----------   --------\n");
  for (int pi = 0; pi < PRIME_COUNT; pi++) {
    double pval = (double)primes[pi];
    double omega = log(pval);
    Cpx dft = dft_at_freq(da, N, omega);
    double cs_err = conj_symmetry_error(da, N, omega);
    printf("%3d  %9.4f   %+.4e   %+.4e   %.4e   %.2e\n", primes[pi], omega,
           dft.re, dft.im, cpx_abs(dft), cs_err);
  }

  printf("\n--- Gram points (reference, not used) ---\n");
  printf(" n   g_n          next zero    spacing\n");
  printf("---  -----------  -----------  -----------\n");
  for (int ni = 0; ni < N; ni++) {
    double g = gram_point(ni);
    printf("%2d   %10.4f    %10.4f    %8.4f\n", ni, g, zeta[ni], zeta[ni] - g);
  }

  return 0;
}
