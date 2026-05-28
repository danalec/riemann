/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   WKB resolvent phase analysis
 * @paper   yamaguchi-rh-2026.tex, §11.12
 * @theorem Theorem III
 * @proof   WKB eigenvectors for free Gram Jacobi
 * @step    4
 *
 * WKB resolvent phase analysis. Paper Appendix D.4, Gap L8.
 * wkb_resolvent.c -- WKB resolvent phase analysis
 *
 * Addresses the L8 conjecture and the k(p) analytic derivation gap.
 *
 * The WKB approximation for the free Gram Jacobi matrix:
 *   u_k(n) ~ sqrt(2/N) sin(n * theta(lambda_k))
 *
 * The resolvent R(z) = (J_free - z)^{-1} has matrix elements:
 *   R_{nm}(z) ~ [psi_< (min(n,m)) * psi_> (max(n,m))] / W
 *
 * where psi_<, psi_> are the decaying/growing WKB solutions
 * and W is the Wronskian.
 *
 * WKB resolvent phase analysis. Paper Appendix D.4, Gap L8.
 * The WKB resolvent phase conjecture (L8):
 *   arg R_{nn}(E + i0) = -pi * S(E)   at energy E
 *
 * This connects the resolvent phase to the scattering phase S(E),
 * which via the explicit formula encodes prime frequencies.
 *
 * This program:
 * 1. Computes WKB eigenvectors for the free Gram Jacobi matrix
 * 2. Computes resolvent matrix elements at prime frequencies
 * 3. Extracts the phase of R_{nn} and tests arg R = -pi*S
 * 4. Derives k(p) from the resolvent phase structure
 *
 * Compile: gcc -Wall -Wextra -Wconversion -Wshadow -Werror -O3
 *          -Isrc -o wkb_resolvent src/wkb_resolvent.c -lm
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

#define MAXN 100

/* Riemann-Siegel theta and derivative */
static double theta_s(double t) {
  if (t <= 1.0)
    return -M_PI / 8.0;
  double u = 1.0 / t;
  double u2 = u * u, u4 = u2 * u2;
  return 0.5 * t * log(t / (2.0 * M_PI)) - 0.5 * t - M_PI / 8.0 + u / 48.0 +
         7.0 * u * u2 / 5760.0 + 31.0 * u4 / 80640.0 +
         127.0 * u4 * u2 / 430080.0;
}

static double theta_p(double t) {
  if (t <= 2.0 * M_PI)
    return 1.0;
  double t2 = t * t;
  return 0.5 * log(t / (2.0 * M_PI)) - 1.0 / (24.0 * t2) +
         7.0 / (960.0 * t2 * t2) + 31.0 / (8064.0 * t2 * t2 * t2);
}

/* WKB eigenvector approximation:
 * u_k(n) ~ sqrt(2/(N * rho(E_k))) * sin((n+1) * theta(E_k))
 * where rho(E) = theta'(E)/pi is the density of states.
 *
 * For the free Gram Jacobi matrix, the WKB approximation gives:
 * u_k(n) ~ C * sin((n+1) * theta(lambda_k) / (N+1))
 */
static void wkb_eigenvector(int N, double lambda, double *u) {
  double phase = theta_s(lambda);
  double rho = theta_p(lambda) / M_PI;
  double norm = sqrt(2.0 / ((double)N * fmax(rho, 0.01)));

  for (int n = 0; n < N; n++) {
    u[n] = norm * sin((n + 1) * phase / (double)(N + 1));
  }

  /* Normalize */
  double sum = 0.0;
  for (int n = 0; n < N; n++)
    sum += u[n] * u[n];
  sum = sqrt(sum);
  if (sum > 1e-15)
    for (int n = 0; n < N; n++)
      u[n] /= sum;
}

/* Build free Gram Jacobi matrix */
static double gram_point(int n) {
  if (n <= 0)
    return 0.0;
  double lo = 10.0, hi = 1e8;
  for (int iter = 0; iter < 80; iter++) {
    double mid = 0.5 * (lo + hi);
    double t = theta_s(mid);
    if (t / M_PI < (double)n)
      lo = mid;
    else
      hi = mid;
  }
  return 0.5 * (lo + hi);
}

static void build_free_jacobi(int N, double *a, double *b) {
  for (int k = 0; k < N; k++) {
    double g_prev = (k > 0) ? gram_point(k - 1) : 1.0;
    double tp = theta_p(fmax(g_prev, 10.0));
    a[k] = g_prev + M_PI / tp;
  }
  for (int k = 0; k < N - 1; k++) {
    double g = gram_point(k);
    double g_next = gram_point(k + 1);
    double tp = theta_p(fmax(g, 10.0));
    b[k] = sqrt(g_next - g) * tp / 2.0;
  }
}

/* Jacobi eigenvalue algorithm (simple cyclic Jacobi) */
static void jacobi_diag(int N, const double *a_in, const double *b_in,
                        double *ev, double *vecs) {
  double J[MAXN * MAXN];
  memset(J, 0, sizeof(J));
  for (int i = 0; i < N; i++)
    J[i * N + i] = a_in[i];
  for (int i = 0; i < N - 1; i++) {
    J[i * N + (i + 1)] = b_in[i];
    J[(i + 1) * N + i] = b_in[i];
  }

  /* Initialize eigenvectors */
  for (int i = 0; i < N; i++) {
    ev[i] = J[i * N + i];
    for (int j = 0; j < N; j++)
      vecs[i * N + j] = (i == j) ? 1.0 : 0.0;
  }

  size_t nn = (size_t)N * (size_t)N;
  double *Jwork = (double *)malloc(nn * sizeof(double));
  memcpy(Jwork, J, nn * sizeof(double));

  int max_iter = 100 * N;
  for (int iter = 0; iter < max_iter; iter++) {
    int p = 0, q = 1;
    double max_val = 0.0;
    for (int i = 0; i < N - 1; i++) {
      for (int j = i + 1; j < N; j++) {
        double v = fabs(Jwork[i * N + j]);
        if (v > max_val) {
          max_val = v;
          p = i;
          q = j;
        }
      }
    }
    if (max_val < 1e-15)
      break;

    double app = Jwork[p * N + p];
    double aqq = Jwork[q * N + q];
    double apq = Jwork[p * N + q];
    double tau = (aqq - app) / (2.0 * apq);
    double t = (tau >= 0) ? 1.0 / (tau + sqrt(1.0 + tau * tau))
                          : 1.0 / (tau - sqrt(1.0 + tau * tau));
    double c = 1.0 / sqrt(1.0 + t * t);
    double s = t * c;

    for (int i = 0; i < N; i++) {
      double dip = Jwork[i * N + p];
      double diq = Jwork[i * N + q];
      Jwork[i * N + p] = c * dip - s * diq;
      Jwork[i * N + q] = s * dip + c * diq;
    }
    for (int i = 0; i < N; i++) {
      double dpi = Jwork[p * N + i];
      double dqi = Jwork[q * N + i];
      Jwork[p * N + i] = c * dpi - s * dqi;
      Jwork[q * N + i] = s * dpi + c * dqi;
    }
    for (int i = 0; i < N; i++) {
      double vip = vecs[i * N + p];
      double viq = vecs[i * N + q];
      vecs[i * N + p] = c * vip - s * viq;
      vecs[i * N + q] = s * vip + c * viq;
    }
  }

  for (int i = 0; i < N; i++)
    ev[i] = Jwork[i * N + i];
  free(Jwork);

  /* Sort by eigenvalue */
  for (int i = 0; i < N - 1; i++) {
    int mi = i;
    for (int j = i + 1; j < N; j++)
      if (ev[j] < ev[mi])
        mi = j;
    if (mi != i) {
      double tmp = ev[i];
      ev[i] = ev[mi];
      ev[mi] = tmp;
      for (int j = 0; j < N; j++) {
        tmp = vecs[j * N + i];
        vecs[j * N + i] = vecs[j * N + mi];
        vecs[j * N + mi] = tmp;
      }
    }
  }
}

/* Compute resolvent diagonal element R_{nn}(z) for z = E + i*eps
 * Using spectral decomposition: R_{nn}(z) = sum_k |u_k(n)|^2 / (lambda_k - z)
 */
static void resolvent_diag(int N, const double *ev, const double *vecs, int n,
                           double E, double eps, double *out_re,
                           double *out_im) {
  double re = 0.0, im = 0.0;
  for (int k = 0; k < N; k++) {
    double u2 = vecs[n * N + k] * vecs[n * N + k];
    double denom = (ev[k] - E) * (ev[k] - E) + eps * eps;
    re += u2 * (ev[k] - E) / denom;
    im += u2 * eps / denom;
  }
  *out_re = re;
  *out_im = im;
}

/* S(E) = E/pi - theta(E)/pi - 1  (Riemann-von Mangoldt) */
static double S_of_E(double E) { return E / M_PI - theta_s(E) / M_PI - 1.0; }

/* Compute k(p) from resolvent phase structure.
 *
 * The key insight: the resolvent R(z) at prime frequency omega = log p
 * has a phase shift delta_p that encodes the perturbation strength.
 *
 * delta_p(E) = arg R_{nn}(E + i0) - arg R_free(E + i0)
 *
 * The perturbation alpha(p) is related to delta_p by:
 * alpha(p) = delta_p / (pi * k(p))
 *
 * So k(p) = delta_p / (pi * alpha(p))
 *
 * From the explicit formula: alpha(p) = -log(p)/(2*pi*sqrt(p))
 * So k(p) = -2*sqrt(p)*delta_p / log(p)
 *
 * If delta_p ~ log(p)/(2*sqrt(p)), then k(p) = 1.
 * The log(p) growth in k(p) comes from the phase accumulation.
 */
static double compute_kp_from_resolvent(int N, const double *ev,
                                        const double *vecs, int p, double E_min,
                                        double E_max, int npts) {
  double omega = log((double)p);
  double dE = (E_max - E_min) / (double)(npts - 1);
  double phase_sum = 0.0;

  for (int i = 0; i < npts; i++) {
    double E = E_min + i * dE;
    double eps = 0.01;

    /* Average over diagonal resolvent elements */
    double re_sum = 0.0, im_sum = 0.0;
    int n_avg = N / 4; /* Average over first quarter */
    for (int n = 0; n < n_avg; n++) {
      double re, im;
      resolvent_diag(N, ev, vecs, n, E, eps, &re, &im);
      re_sum += re;
      im_sum += im;
    }
    re_sum /= (double)n_avg;
    im_sum /= (double)n_avg;

    /* Phase of resolvent */
    double phase = atan2(im_sum, re_sum);

    /* Weighted by heat kernel */
    double weight = exp(-E / 100.0);
    phase_sum += weight * phase * sin(omega * E);
  }

  phase_sum *= dE;

  /* k(p) = |phase_sum| * normalization */
  double kp = fabs(phase_sum) * 2.0 * sqrt((double)p) / omega;
  return kp;
}

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

  printf("================================================================\n");
  printf("  WKB Resolvent Phase Analysis\n");
  printf("  L8 Conjecture and k(p) Analytic Derivation\n");
  printf(
      "================================================================\n\n");

  int primes[] = {5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59};
  int np = sizeof(primes) / sizeof(primes[0]);

  int N = 30;
  double a[MAXN], b[MAXN];
  double ev[MAXN], vecs[MAXN * MAXN];
  double wkb_u[MAXN];

  printf("TEST 1: WKB eigenvector vs exact eigenvector\n\n");

  build_free_jacobi(N, a, b);
  jacobi_diag(N, a, b, ev, vecs);

  printf("  %4s  %10s  %10s  %10s  %10s  %10s\n", "k", "lambda_k",
         "||u_exact||", "||u_wkb||", "<u_e|u_w>", "angle_deg");
  printf("  %4s  %10s  %10s  %10s  %10s  %10s\n", "---", "--------",
         "---------", "--------", "--------", "--------");

  for (int k = 0; k < N && k < 10; k++) {
    wkb_eigenvector(N, ev[k], wkb_u);

    double dot = 0.0;
    double norm_wkb = 0.0;
    for (int n = 0; n < N; n++) {
      dot += vecs[n * N + k] * wkb_u[n];
      norm_wkb += wkb_u[n] * wkb_u[n];
    }
    norm_wkb = sqrt(norm_wkb);
    double angle =
        fabs(dot) > 1e-15 ? acos(fabs(dot) / norm_wkb) * 180.0 / M_PI : 90.0;

    printf("  %4d  %10.4f  %10.6f  %10.6f  %10.6f  %10.2f\n", k, ev[k], 1.0,
           norm_wkb, dot, angle);
  }

  printf("\n\nTEST 2: Resolvent phase vs S(E)\n\n");

  printf("  Testing arg R_{nn}(E+i0) = -pi*S(E)\n\n");

  printf("  %8s  %10s  %10s  %10s  %10s  %8s\n", "E", "arg R", "-pi*S(E)",
         "Re R", "Im R", "error");
  printf("  %8s  %10s  %10s  %10s  %10s  %8s\n", "------", "------", "------",
         "------", "------", "------");

  double E_test[] = {20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0};
  int nE = sizeof(E_test) / sizeof(E_test[0]);
  int n_mid = N / 2;

  for (int ei = 0; ei < nE; ei++) {
    double E = E_test[ei];
    double re, im;
    resolvent_diag(N, ev, vecs, n_mid, E, 0.01, &re, &im);
    double arg_R = atan2(im, re);
    double neg_pi_S = -M_PI * S_of_E(E);
    double error = fabs(arg_R - neg_pi_S);

    /* Wrap error to [-pi, pi] */
    while (error > M_PI)
      error -= 2.0 * M_PI;
    error = fabs(error);

    printf("  %8.1f  %10.4f  %10.4f  %10.6f  %10.6f  %8.4f\n", E, arg_R,
           neg_pi_S, re, im, error);
  }

  printf("\n\nTEST 3: k(p) from resolvent phase\n\n");

  printf("  %4s  %8s  %10s  %10s  %8s\n", "p", "log(p)", "k(p)_WKB", "k(p)_emp",
         "ratio");
  printf("  %4s  %8s  %10s  %10s  %8s\n", "---", "------", "--------",
         "--------", "------");

  double k_wkb[MAXN];
  double x_wkb[MAXN];
  int nk = 0;

  for (int i = 0; i < np; i++) {
    int p = primes[i];
    double log_p = log((double)p);
    double kp_wkb =
        compute_kp_from_resolvent(N, ev, vecs, p, 15.0, 135.0, 5000);
    double kp_emp = 0.27 + 0.66 * log_p;
    double ratio = (fabs(kp_emp) > 1e-15) ? kp_wkb / kp_emp : 0.0;

    k_wkb[nk] = kp_wkb;
    x_wkb[nk] = log_p;
    nk++;

    printf("  %4d  %8.3f  %10.6f  %10.6f  %8.4f\n", p, log_p, kp_wkb, kp_emp,
           ratio);
  }

  printf("\n\nTEST 4: Fit k(p) = A + B log p to WKB results\n\n");

  double sum_x = 0.0, sum_y = 0.0, sum_xx = 0.0, sum_xy = 0.0;
  for (int i = 0; i < nk; i++) {
    sum_x += x_wkb[i];
    sum_y += k_wkb[i];
    sum_xx += x_wkb[i] * x_wkb[i];
    sum_xy += x_wkb[i] * k_wkb[i];
  }
  double mean_x = sum_x / (double)nk;
  double mean_y = sum_y / (double)nk;
  double B = (sum_xy - (double)nk * mean_x * mean_y) /
             (sum_xx - (double)nk * mean_x * mean_x);
  double A = mean_y - B * mean_x;

  double ss_res = 0.0, ss_tot = 0.0;
  for (int i = 0; i < nk; i++) {
    double pred = A + B * x_wkb[i];
    double err = k_wkb[i] - pred;
    ss_res += err * err;
    ss_tot += (k_wkb[i] - mean_y) * (k_wkb[i] - mean_y);
  }
  double R2 = (ss_tot > 1e-15) ? 1.0 - ss_res / ss_tot : 0.0;

  printf("  WKB fitted:  k(p) = %.4f + %.4f log p  (R^2 = %.6f)\n", A, B, R2);
  printf("  Empirical:   k(p) = 0.2700 + 0.6600 log p\n\n");

  printf("  A: WKB = %.4f, empirical = 0.27, diff = %+.4f\n", A, A - 0.27);
  printf("  B: WKB = %.4f, empirical = 0.66, diff = %+.4f\n", B, B - 0.66);

  printf(
      "\n================================================================\n");
  printf("  CONCLUSION\n");
  printf(
      "================================================================\n\n");

  printf("  1. WKB eigenvector approximation:\n");
  printf("     Accuracy: see TEST 1 (angle between exact and WKB)\n\n");
  printf("  2. Resolvent phase arg R vs -pi*S(E):\n");
  printf("     See TEST 2 (error column)\n\n");
  // WKB resolvent phase analysis. Paper Appendix D.4, Gap L8.
  printf("  3. k(p) from WKB resolvent:\n");
  printf("     k(p) = %.4f + %.4f log p  (R^2 = %.6f)\n\n", A, B, R2);

  if (R2 > 0.9) {
    // WKB resolvent phase analysis. Paper Appendix D.4, Gap L8.
    printf("  PASS: WKB resolvent confirms k(p) = A + B log p.\n");
    if (fabs(A - 0.27) < 0.5 && fabs(B - 0.66) < 0.5) {
      printf("  Coefficients are within 0.5 of empirical values.\n");
    }
  } else {
    // WKB resolvent phase analysis. Paper Appendix D.4, Gap L8.
    printf("  The WKB resolvent phase model needs refinement.\n");
    printf("  Possible improvements:\n");
    printf("  - Use exact eigenvectors instead of WKB approx\n");
    printf("  - Include off-diagonal resolvent elements\n");
    printf("  - Account for boundary effects at small N\n");
  }

  return 0;
}
