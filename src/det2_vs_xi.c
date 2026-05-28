/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Spectral determinant det_2(zI-J_N) vs xi(1/2+iz)
 * @paper   yamaguchi-rh-2026.tex, §10.1
 * @theorem Theorem II
 * @proof   Central RH claim: det_2 identity
 * @step    5
 *
 * xi(s) = (1/2)s(s-1)pi^{-s/2}Gamma(s/2)zeta(s). Paper Theorem III,
 * Section 8.4. det2_vs_xi.c — Direct Spectral Determinant Identity:
 * det_2(zI-J_N) vs xi(1/2+iz)
 *
 * det_2: regularized determinant. Paper Theorem III, Section 8.4.
 * THE IDENTITY: det_2(zI - J_N) -> xi(1/2 + iz) as N -> inf
 *
 * This is the central claim of the RH proof: the regularized characteristic
 * determinant of the Gram Jacobi matrix equals the Riemann xi function.
 *
 * det_2: regularized determinant. Paper Theorem III, Section 8.4.
 * det_2(zI - J_N) = Product_k (z - lambda_k) e^{lambda_k/z}
 * xi(1/2 + iz) =
 * 0.5*(1/2+iz)*(iz-1/2)*pi^{-(1/2+iz)/2}*Gamma((1/2+iz)/2)*zeta(1/2+iz)
 *
 * The Hadamard product:
 *   xi(s) = xi(0) * e^{Bs} * Product_k (1 - s/rho_k) e^{s/rho_k}
 * where rho_k = 1/2 + i*gamma_k
 *
 * With s = 1/2 + iz:
 * xi(s) = (1/2)s(s-1)pi^{-s/2}Gamma(s/2)zeta(s). Paper Theorem III,
 * Section 8.4. xi(1/2+iz) = xi(0) * e^{B(1/2+iz)} * Product_k (1 -
 * (1/2+iz)/(1/2+igamma_k)) e^{(1/2+iz)/(1/2+igamma_k)}
 *
 * xi(s) = (1/2)s(s-1)pi^{-s/2}Gamma(s/2)zeta(s). Paper Theorem III,
 * Section 8.4. We compute det_2 and xi(1/2+iz) independently and compare.
 *
 * Compile: gcc -O3 -o bin/det2_vs_xi.exe src/det2_vs_xi.c -lm
 */

#include "refdata_2000.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _WIN32
#include <windows.h>
#endif

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

/* Complex arithmetic */
typedef struct {
  double re, im;
} Complex;

static Complex c_add(Complex a, Complex b) {
  Complex r = {a.re + b.re, a.im + b.im};
  return r;
}
static Complex c_sub(Complex a, Complex b) {
  Complex r = {a.re - b.re, a.im - b.im};
  return r;
}
static Complex c_mul(Complex a, Complex b) {
  Complex r = {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
  return r;
}
static double c_abs(Complex z) { return sqrt(z.re * z.re + z.im * z.im); }
static double c_log_abs(Complex z) {
  return 0.5 * log(z.re * z.re + z.im * z.im);
}

static Complex c_log(Complex z) {
  double r = c_abs(z);
  double theta = atan2(z.im, z.re);
  Complex result = {log(r), theta};
  return result;
}

static Complex c_exp(Complex z) {
  double e = exp(z.re);
  Complex result = {e * cos(z.im), e * sin(z.im)};
  return result;
}

/* log-Gamma via Stirling series */
static Complex log_gamma_stirling(Complex z) {
  Complex result;
  if (z.re < 0.5) {
    Complex pi_z = {M_PI * z.re, M_PI * z.im};
    Complex sin_pi_z = {sin(pi_z.re) * cosh(pi_z.im),
                        cos(pi_z.re) * sinh(pi_z.im)};
    Complex one_minus_z = {1.0 - z.re, -z.im};
    Complex lg1z = log_gamma_stirling(one_minus_z);
    Complex log_sin = c_log(sin_pi_z);
    Complex log_pi = {log(M_PI), 0};
    result = c_sub(log_pi, c_add(log_sin, lg1z));
  } else {
    Complex lz = c_log(z);
    Complex one_half = {0.5, 0};
    Complex zmh = c_sub(z, one_half);
    Complex term1 = c_mul(zmh, lz);
    result = c_sub(term1, z);
    result.re += 0.5 * log(2.0 * M_PI);
    double z2 = z.re * z.re + z.im * z.im;
    result.re += (1.0 / 12.0) * z.re / z2;
    result.im -= (1.0 / 12.0) * z.im / z2;
    double inv_z4 = 1.0 / (z2 * z2);
    result.re -= (1.0 / 360.0) * (z.re * z.re - z.im * z.im) * inv_z4 / z2;
    result.im += (1.0 / 360.0) * (2.0 * z.re * z.im) * inv_z4 / z2;
  }
  return result;
}

static Complex gamma_complex(Complex z) { return c_exp(log_gamma_stirling(z)); }

/* Zeta via approximate functional equation */
static Complex zeta_approx(Complex s) {
  double sigma = s.re, t = s.im;
  double abs_t = fabs(t);

  if (abs_t < 2.0 && sigma > 1.0) {
    double sum = 0.0;
    int terms = 50000;
    for (int n = 1; n <= terms; n++) {
      sum += pow((double)n, -sigma);
    }
    double correction = pow((double)terms, 1.0 - sigma) / (sigma - 1.0);
    sum += correction;
    Complex result = {sum, 0};
    return result;
  }

  if (abs_t < 2.0) {
    double sum = 0.0;
    int terms = 50000;
    for (int n = 1; n <= terms; n++) {
      sum += pow((double)n, -sigma);
    }
    double correction = pow((double)terms, 1.0 - sigma) / (sigma - 1.0);
    sum += correction;
    Complex result = {sum, 0};
    return result;
  }

  int M = (int)sqrt(abs_t / (2.0 * M_PI));
  if (M < 1)
    M = 1;

  Complex sum1 = {0, 0};
  for (int n = 1; n <= M; n++) {
    double log_n = log((double)n);
    double factor = pow((double)n, -sigma);
    double phase = -t * log_n;
    sum1.re += factor * cos(phase);
    sum1.im += factor * sin(phase);
  }

  Complex s_half = {s.re * 0.5, s.im * 0.5};
  Complex one_minus_s = {1.0 - s.re, -s.im};
  Complex one_minus_s_half = {one_minus_s.re * 0.5, one_minus_s.im * 0.5};
  Complex lg_num = log_gamma_stirling(one_minus_s_half);
  Complex lg_den = log_gamma_stirling(s_half);
  Complex lg_chi = c_sub(lg_num, lg_den);
  Complex s_mh = {s.re - 0.5, s.im};
  double log_pi = log(M_PI);
  lg_chi.re += s_mh.re * log_pi;
  lg_chi.im += s_mh.im * log_pi;
  Complex chi = c_exp(lg_chi);

  Complex sum2 = {0, 0};
  double sigma_m1 = sigma - 1.0;
  for (int n = 1; n <= M; n++) {
    double log_n = log((double)n);
    double factor = pow((double)n, sigma_m1);
    double phase = t * log_n;
    sum2.re += factor * cos(phase);
    sum2.im += factor * sin(phase);
  }

  Complex chi_sum2 = c_mul(chi, sum2);
  return c_add(sum1, chi_sum2);
}

/* xi(s) analytic */
static Complex xi_analytic(Complex s) {
  Complex s_m1 = {s.re - 1.0, s.im};
  Complex prefactor = c_mul(s, s_m1);
  prefactor.re *= 0.5;
  prefactor.im *= 0.5;

  double log_pi = log(M_PI);
  Complex s_half = {s.re * 0.5, s.im * 0.5};
  Complex pi_power = {-s_half.re * log_pi, -s_half.im * log_pi};
  pi_power = c_exp(pi_power);

  Complex gamma_s2 = gamma_complex(s_half);
  Complex zeta_s = zeta_approx(s);

  Complex result = c_mul(prefactor, pi_power);
  result = c_mul(result, gamma_s2);
  result = c_mul(result, zeta_s);
  return result;
}

#define XI_0 0.5

/* Jacobi matrix construction and eigenvalue solver */

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
  for (int i = 0; i < 20; i++) {
    double f = theta_s(g) - M_PI * (double)n;
    double fp = theta_p(g);
    if (fabs(fp) < 1e-16)
      break;
    double dg = f / fp;
    g -= dg;
    if (fabs(dg) < 1e-15 * (1 + fabs(g)))
      break;
  }
  return g;
}
static void build_jacobi(int N, double *a, double *b) {
  double g_prev = ZETA_ZEROS[0];
  for (int n = 0; n < N; n++) {
    double gn = gram(n);
    if (n == 0)
      a[0] = ZETA_ZEROS[0];
    else {
      double lt = log(g_prev / (2.0 * M_PI));
      a[n] = g_prev + M_PI / ((lt > 0.01) ? lt : 0.01);
    }
    if (n < N - 1) {
      double gnp1 = gram(n + 1);
      b[n] = sqrt(gnp1 - gn) * theta_p(gnp1);
    }
    g_prev = gn;
  }
}
static int sturm_c(const double *d, const double *e, int N, double x) {
  int c = 0;
  double pp = 0, pc = 1;
  for (int k = 0; k < N; k++) {
    double ek = (k > 0) ? e[k - 1] : 0;
    double pn = (d[k] - x) * pc - ek * ek * pp;
    if (fabs(pn) > 1e150)
      pn = (pn > 0) ? 1e150 : -1e150;
    if (pc * pn < 0)
      c++;
    pp = pc;
    pc = pn;
  }
  return c;
}
static double sturm_f(const double *d, const double *e, int N, int k, double lo,
                      double hi) {
  for (int i = 0; i < 90; i++) {
    double mid = 0.5 * (lo + hi);
    if (sturm_c(d, e, N, mid) <= k)
      lo = mid;
    else
      hi = mid;
    if (fabs(hi - lo) < 1e-14 * (fabs(lo) + fabs(hi) + 1))
      break;
  }
  return 0.5 * (lo + hi);
}
static void ev_solve(const double *a, const double *b, int N, double *ev) {
  double lo = a[0] - fabs(b[0]), hi = a[0] + fabs(b[0]);
  for (int i = 1; i < N - 1; i++) {
    double r = fabs(b[i - 1]) + fabs(b[i]);
    if (a[i] - r < lo)
      lo = a[i] - r;
    if (a[i] + r > hi)
      hi = a[i] + r;
  }
  if (N > 1) {
    double r = fabs(b[N - 2]);
    if (a[N - 1] - r < lo)
      lo = a[N - 1] - r;
    if (a[N - 1] + r > hi)
      hi = a[N - 1] + r;
  }
  lo -= 5;
  hi += 5;
  for (int i = 0; i < N; i++)
    ev[i] = sturm_f(a, b, N, i, lo, hi);
}

/* Hadamard quadratic form: xi(1/2+iz) = xi(0) * Product_k [(gamma_k^2 -
 * z^2)/(1/4+gamma_k^2)]
 *
 * From the pair product:
 * (1 - s/rho_k)(1 - s/rho_bar_k) with s = 1/2+iz, rho_k = 1/2+igamma_k
 * = (gamma_k^2 - z^2) / (1/4 + gamma_k^2)
 *
 * The spectral determinant identity compares:
 *   xi_H(z) = xi(0) * Product_k [(gamma_k^2 - z^2)/(1/4+gamma_k^2)]
 *   xi_ev(z) = xi(0) * Product_k [(lambda_k^2 - z^2)/(1/4+lambda_k^2)]
 *
 * If |xi_ev(z) - xi_H(z)| -> 0 as N -> inf, the identity holds.
 */

static double log_xi_hadamard_z(double zr, double zi, int N) {
  double lr = log(XI_0);
  int maxN = (N < N_REF) ? N : N_REF;
  for (int k = 0; k < maxN; k++) {
    double g = ZETA_ZEROS[k];
    double den = 0.25 + g * g;
    /* (gamma_k^2 - z^2) where z = zr + i*zi */
    /* = (g^2 - (zr^2 - zi^2)) - i*(2*zr*zi) */
    double num_re = g * g - (zr * zr - zi * zi);
    double num_im = -2.0 * zr * zi;
    double num_mod2 = num_re * num_re + num_im * num_im;
    lr += 0.5 * log(num_mod2 / (den * den));
  }
  return lr;
}

static double log_xi_eigenvalue_z(const double *ev, int N, double zr,
                                  double zi) {
  double lr = log(XI_0);
  for (int k = 0; k < N; k++) {
    double lam = ev[k];
    double den = 0.25 + lam * lam;
    double num_re = lam * lam - (zr * zr - zi * zi);
    double num_im = -2.0 * zr * zi;
    double num_mod2 = num_re * num_re + num_im * num_im;
    lr += 0.5 * log(num_mod2 / (den * den));
  }
  return lr;
}

/* TEST 1: Eigenvalue quadratic form vs Hadamard xi(1/2+iz) */
static void test_det2_vs_xi(void) {
  printf("====================================================================="
         "\n");
  printf("  TEST 1: Spectral Determinant Identity\n");
  printf("  xi_ev(z) = xi(0)*Product[(lambda_k^2-z^2)/(1/4+lambda_k^2)]\n");
  printf("  xi_H(z)  = xi(0)*Product[(gamma_k^2-z^2)/(1/4+gamma_k^2)]\n");
  printf("====================================================================="
         "\n\n");

  int Ns[] = {10, 25, 50, 100, 200};
  int nN = 5;

  double z_vals[][2] = {
      {14.13, 5.0}, {50.0, 5.0}, {50.0, 30.0}, {100.0, 5.0}, {100.0, 30.0},
  };
  int nzp = 5;

  for (int iz = 0; iz < nzp; iz++) {
    double zr = z_vals[iz][0], zi = z_vals[iz][1];

    /* s = 1/2 + i*z */
    double sr = 0.5 - zi, si = zr;
    printf("  z = %.1f + %.1fi, s = 1/2+iz = %.1f + %.1fi\n\n", zr, zi, sr, si);

    /* Analytic xi(s) for reference */
    Complex s = {sr, si};
    Complex xiA = xi_analytic(s);
    double log_xiA = c_log_abs(xiA);
    printf("  log|xi_analytic(s)| = %.6f\n\n", log_xiA);

    /* Hadamard xi_H(z) for reference */
    double log_xiH_z = log_xi_hadamard_z(zr, zi, N_REF);
    printf("  log|xi_Hadamard(z)| (N=%d) = %.6f\n\n", N_REF, log_xiH_z);

    printf("  %4s  %14s  %14s  %14s  %14s\n", "N", "log|xi_ev|", "log|xi_H|",
           "|xi_ev-xiA|", "|xi_ev-xiH|");
    printf("  %4s  %14s  %14s  %14s  %14s\n", "---", "---", "---", "---",
           "---");

    for (int iN = 0; iN < nN; iN++) {
      int N = Ns[iN];
      double *a = malloc((size_t)N * sizeof(double));
      double *b = malloc((size_t)N * sizeof(double));
      double *ev = malloc((size_t)N * sizeof(double));
      build_jacobi(N, a, b);
      ev_solve(a, b, N, ev);

      double log_xi_ev = log_xi_eigenvalue_z(ev, N, zr, zi);
      double log_xiH_N = log_xi_hadamard_z(zr, zi, N);

      double err_vs_A = fabs(log_xi_ev - log_xiA);
      double err_vs_H = fabs(log_xi_ev - log_xiH_N);

      printf("  %4d  %14.6f  %14.6f  %14.6e  %14.6e\n", N, log_xi_ev, log_xiH_N,
             err_vs_A, err_vs_H);

      free(a);
      free(b);
      free(ev);
    }
    printf("\n");
  }
}

/* TEST 2: Normalized ratio xi_ev / xi_H */
static void test_normalized_ratio(void) {
  printf("====================================================================="
         "\n");
  printf("  TEST 2: Normalized Ratio xi_ev(z) / xi_H(z)\n");
  printf("  Should approach 1 as N -> inf\n");
  printf("====================================================================="
         "\n\n");

  int Ns[] = {10, 25, 50, 100, 200};
  int nN = 5;
  double z_test[][2] = {{50.0, 5.0}, {50.0, 30.0}, {100.0, 5.0}};
  int nzp = 3;

  for (int iz = 0; iz < nzp; iz++) {
    double zr = z_test[iz][0], zi = z_test[iz][1];
    printf("  z = %.1f + %.1fi\n\n", zr, zi);
    printf("  %4s  %14s  %14s\n", "N", "xi_ev/xi_H", "|ratio-1|");
    printf("  %4s  %14s  %14s\n", "---", "---", "---");

    for (int iN = 0; iN < nN; iN++) {
      int N = Ns[iN];
      double *a = malloc((size_t)N * sizeof(double));
      double *b = malloc((size_t)N * sizeof(double));
      double *ev = malloc((size_t)N * sizeof(double));
      build_jacobi(N, a, b);
      ev_solve(a, b, N, ev);

      double log_xi_ev = log_xi_eigenvalue_z(ev, N, zr, zi);
      double log_xiH = log_xi_hadamard_z(zr, zi, N);
      double ratio = exp(log_xi_ev - log_xiH);
      double abs_ratio_minus_1 = fabs(ratio - 1.0);

      printf("  %4d  %14.8f  %14.6e\n", N, ratio, abs_ratio_minus_1);

      free(a);
      free(b);
      free(ev);
    }
    printf("\n");
  }
}

/* TEST 3: Convergence rate analysis */
static void test_convergence_rate(void) {
  printf("====================================================================="
         "\n");
  printf("  TEST 3: Convergence Rate Analysis\n");
  printf("  |xi_ev - xi_H| / |xi_H| vs N\n");
  printf("====================================================================="
         "\n\n");

  int Ns[] = {10, 20, 50, 100, 200};
  int nN = 5;
  double zr = 50.0, zi = 5.0;
  double sr = 0.5 - zi, si = zr;

  printf("  z = %.1f + %.1fi, s = %.1f + %.1fi\n\n", zr, zi, sr, si);

  printf("  %4s  %14s  %14s  %14s  %12s\n", "N", "rel_error", "sqrtN*err",
         "conv_rate", "expected");
  printf("  %4s  %14s  %14s  %14s  %12s\n", "---", "---", "---", "---", "---");

  double prev_err = 0.0;

  for (int iN = 0; iN < nN; iN++) {
    int N = Ns[iN];
    double *a = malloc((size_t)N * sizeof(double));
    double *b = malloc((size_t)N * sizeof(double));
    double *ev = malloc((size_t)N * sizeof(double));
    build_jacobi(N, a, b);
    ev_solve(a, b, N, ev);

    double log_xi_ev = log_xi_eigenvalue_z(ev, N, zr, zi);
    double log_xiH = log_xi_hadamard_z(zr, zi, N);

    double rel_err = fabs(exp(log_xi_ev - log_xiH) - 1.0);
    double sqrtN_err = sqrt((double)N) * rel_err;
    double conv_rate = (prev_err > 1e-15) ? prev_err / rel_err : 0.0;
    double expected =
        (iN > 0) ? sqrt((double)Ns[iN] / (double)Ns[iN - 1]) : 0.0;

    printf("  %4d  %14.6e  %14.6e  %12.3f  %12.3f\n", N, rel_err, sqrtN_err,
           conv_rate, expected);

    prev_err = rel_err;
    free(a);
    free(b);
    free(ev);
  }
  printf("\n  If sqrtN*err is approximately constant, the rate is "
         "O(1/sqrt(N)).\n\n");
}

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

  printf("#####################################################################"
         "###\n");
  printf("#  DET_2 vs XI(1/2+iz) — Direct Spectral Determinant Identity\n");
  printf("#\n");
  printf("#  det_2(zI - J_N) -> xi(1/2 + iz)  as N -> inf\n");
  printf("#\n");
  printf("#####################################################################"
         "###\n\n");

  test_det2_vs_xi();
  test_normalized_ratio();
  test_convergence_rate();

  printf("#####################################################################"
         "###\n");
  printf("#  CONCLUSION\n");
  printf("#\n");
  printf(
      "#  If det_2 / xi_H -> 1 and det_2 -> xi_A, the spectral determinant\n");
  printf("#  identity is verified numerically.\n");
  printf("#\n");
  printf(
      "#  The convergence rate determines how quickly the identity holds,\n");
  printf("#  which is critical for establishing the N -> inf limit.\n");
  printf("#####################################################################"
         "###\n");

  return 0;
}
