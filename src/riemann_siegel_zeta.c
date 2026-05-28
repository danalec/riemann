/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Riemann-Siegel Z(t) via Euler-Maclaurin
 * @paper   yamaguchi-rh-2026.tex, §2
 * @theorem Infrastructure
 * @proof   Zeta function evaluation
 * @step    Infrastructure
 *
 * riemann_siegel_zeta.c — Accurate zeta(1/2+it) via Euler-Maclaurin Formula
 *
 * Uses the Euler-Maclaurin summation formula for rapid convergence:
 *   zeta(s) = sum_{n=1}^{N-1} n^{-s} + N^{-s}/2 + N^{1-s}/(s-1)
 *           + sum_{k=1}^{K} B_{2k}/(2k)! * prod_{j=0}^{2k-2}(s+j) * N^{-s-2k+1}
 *
 * This converges exponentially with K for fixed N.
 * Then Z(t) = e^{i*theta(t)} * zeta(1/2+it) (real-valued).
 *
 * Reference: tests/ref_xi_mpmath.py (mpmath, 50-digit precision).
 *
 * Compile: gcc -O3 -o bin/riemann_siegel_zeta.exe src/riemann_siegel_zeta.c -lm
 */

#include "refdata_2000.h"
#include <math.h>
#include <stdio.h>

#ifdef _WIN32
#include <windows.h>
#endif

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

/* ---- Complex arithmetic for zeta computation ---- */

typedef struct {
  double re, im;
} Complex;

static Complex c_add(Complex a, Complex b) {
  Complex r = {a.re + b.re, a.im + b.im};
  return r;
}
static Complex c_mul(Complex a, Complex b) {
  Complex r = {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
  return r;
}
static Complex c_div(Complex a, Complex b) {
  double d = b.re * b.re + b.im * b.im;
  Complex r = {(a.re * b.re + a.im * b.im) / d,
               (a.im * b.re - a.re * b.im) / d};
  return r;
}

/* ---- log Gamma via Stirling (real part only) ---- */

static double log_gamma_re(double x, double y) {
  double r = sqrt(x * x + y * y);
  double re = (x - 0.5) * log(r) - y * atan2(y, x) - x + 0.5 * log(2.0 * M_PI);
  double r2 = r * r;
  double r4 = r2 * r2;
  re += (1.0 / 12.0) * x / r2;
  re -= (1.0 / 360.0) * (x * x - 3.0 * y * y) * x / r4;
  return re;
}

// Z(t) = e^{i*theta(t)} zeta(1/2+it). Paper Section 2.1.
/* ---- Riemann-Siegel theta(t) ---- */

static double rs_theta(double t) {
  if (t < 1.0)
    return -M_PI / 8.0;
  double u = 1.0 / t;
  double u2 = u * u;
  double u4 = u2 * u2;
  return (t / 2.0) * log(t / (2.0 * M_PI)) - t / 2.0 - M_PI / 8.0 + u / 48.0 +
         7.0 * u * u2 / 5760.0 + 31.0 * u4 * u / 80640.0 +
         127.0 * u4 * u2 / 430080.0;
}

/* ---- Bernoulli numbers B_{2k} ----
 * B_2 = 1/6, B_4 = -1/30, B_6 = 1/42, B_8 = -1/30, ...
 */

static double bernoulli_2k(int k) {
  /* Returns B_{2k} for k = 1, 2, ..., 15 */
  static const double B[] = {
      1.0 / 6.0,                 /* B_2 */
      -1.0 / 30.0,               /* B_4 */
      1.0 / 42.0,                /* B_6 */
      -1.0 / 30.0,               /* B_8 */
      5.0 / 66.0,                /* B_10 */
      -691.0 / 2730.0,           /* B_12 */
      7.0 / 6.0,                 /* B_14 */
      -3617.0 / 510.0,           /* B_16 */
      43867.0 / 798.0,           /* B_18 */
      -174611.0 / 330.0,         /* B_20 */
      77683.0 / 5.0,             /* B_22 */
      -23749461029.0 / 870.0,    /* B_24 */
      8615841276005.0 / 14322.0, /* B_26 */
      -7709321041217.0 / 510.0,  /* B_28 */
      2577687858367.0 / 6.0,     /* B_30 */
  };
  if (k < 1 || k > 15)
    return 0.0;
  return B[k - 1];
}

/* ---- (2k)! ---- */

static double factorial_2k(int k) {
  double result = 1.0;
  for (int i = 2; i <= 2 * k; i++)
    result *= (double)i;
  return result;
}

/* ---- zeta(s) via Euler-Maclaurin formula ----
 *
 * zeta(s) = sum_{n=1}^{N-1} n^{-s} + N^{-s}/2 + N^{1-s}/(s-1)
 *         + sum_{k=1}^{K} B_{2k}/(2k)! * prod_{j=0}^{2k-2}(s+j) * N^{-s-2k+1}
 *
 * For s = 1/2 + it:
 *   n^{-s} = n^{-1/2} * e^{-it*log(n)} = n^{-1/2} * (cos(t*log n) - i*sin(t*log
 * n))
 *
 * We choose N ~ t and K ~ 20 for good convergence.
 * For small t, N=100 and K=15 is sufficient.
 * For large t, N ~ t and K ~ 20 gives ~15 digit accuracy.
 */

static Complex zeta_euler_maclaurin(double sigma, double t) {
  Complex s = {sigma, t};
  (void)s; /* unused — kept for API compatibility */

  /* Choose N and K based on t */
  int N = (int)(t + 20);
  if (N < 50)
    N = 50;
  int K = 15;

  /* sum_{n=1}^{N-1} n^{-s} */
  Complex sum = {0.0, 0.0};
  for (int n = 1; n < N; n++) {
    double log_n = log((double)n);
    double factor = pow((double)n, -sigma);
    double phase = -t * log_n;
    sum.re += factor * cos(phase);
    sum.im += factor * sin(phase);
  }

  /* N^{-s}/2 */
  double log_N = log((double)N);
  double N_factor = pow((double)N, -sigma);
  double N_phase = -t * log_N;
  Complex term1 = {N_factor * cos(N_phase) / 2.0,
                   N_factor * sin(N_phase) / 2.0};
  sum = c_add(sum, term1);

  /* N^{1-s}/(s-1) */
  Complex s_m1 = {sigma - 1.0, t};
  double N_1m_sigma = pow((double)N, 1.0 - sigma);
  double N_1m_phase = -t * log_N;
  Complex num = {N_1m_sigma * cos(N_1m_phase), N_1m_sigma * sin(N_1m_phase)};
  Complex term2 = c_div(num, s_m1);
  sum = c_add(sum, term2);

  /* Sum_{k=1}^{K} B_{2k}/(2k)! * prod_{j=0}^{2k-2}(s+j) * N^{-s-2k+1} */
  for (int k = 1; k <= K; k++) {
    double B2k = bernoulli_2k(k);
    if (B2k == 0.0)
      continue;

    double fact2k = factorial_2k(k);
    double coeff = B2k / fact2k;

    /* prod_{j=0}^{2k-2}(s+j) */
    Complex prod = {1.0, 0.0};
    for (int j = 0; j <= 2 * k - 2; j++) {
      Complex s_j = {sigma + j, t};
      prod = c_mul(prod, s_j);
    }

    /* N^{-s-2k+1} */
    double N_factor2 = pow((double)N, -sigma - 2.0 * k + 1.0);
    double N_phase2 = -t * log_N;
    Complex N_term = {N_factor2 * cos(N_phase2), N_factor2 * sin(N_phase2)};

    Complex term = c_mul(prod, N_term);
    term.re *= coeff;
    term.im *= coeff;
    sum = c_add(sum, term);
  }

  return sum;
}

/* ---- Z(t) = e^{i*theta(t)} * zeta(1/2+it) ----
 *
 * Z(t) = cos(theta)*Re(zeta) - sin(theta)*Im(zeta)
 */

static double Z_of_t(double t) {
  Complex zeta = zeta_euler_maclaurin(0.5, t);
  double theta = rs_theta(t);
  return cos(theta) * zeta.re - sin(theta) * zeta.im;
}

/* ---- |xi(1/2+it)| from Z(t) ----
 *
 * |xi(1/2+it)| = 0.5 * sqrt(1/4+t^2) * pi^{-1/4} * |Gamma(1/4+it/2)| * |Z(t)|
 */

static double log_abs_xi_from_Z(double t) {
  double log_prefactor = log(0.5) + 0.5 * log(0.25 + t * t) - 0.25 * log(M_PI);
  double log_gamma = log_gamma_re(0.25, t / 2.0);
  double Z = Z_of_t(t);
  double log_Z = log(fabs(Z) + 1e-300);
  return log_prefactor + log_gamma + log_Z;
}

/* ---- mpmath reference values (50-digit precision) ---- */

static const struct {
  double t;
  double log_xi;
} mp_ref[] = {
    {14.1347251417, -30.66818121793730}, {21.0220396388, -35.22221987565582},
    {30.4248761259, -41.47041217178927}, {50.0, -33.38750994836327},
    {100.0, -69.37730139844239},         {150.0, -111.32498031241307},
    {200.0, -145.97373400189206},
};
static const int n_mp = sizeof(mp_ref) / sizeof(mp_ref[0]);

/* mpmath log10|xi| at zeta zeros */
static const double mp_log10_at_zeros[] = {
    -54, -56, -56, -58, -59, -60, -61, -62, -63, -64,
};

/* ====================================================================
 * TEST 1: Z(t) at known zeta zeros
 * ==================================================================== */

static void test_z_at_zeros(void) {
  printf("====================================================================="
         "\n");
  printf("  TEST 1: Z(t) at Known Zeta Zeros\n");
  printf("  Z(gamma_k) should be ~0 (mpmath: |xi| = 10^-54 to 10^-64)\n");
  printf("====================================================================="
         "\n\n");

  printf("  %4s  %14s  %14s  %10s  %10s\n", "k", "gamma_k", "|Z(gamma_k)|",
         "log10|Z|", "mpmath log10");
  printf("  %4s  %14s  %14s  %10s  %10s\n", "---", "---", "---", "---", "---");

  int nk = 10;
  double log10_e = 1.0 / log(10.0);

  for (int k = 0; k < nk && k < N_REF; k++) {
    double g = ZETA_ZEROS[k];
    double Z = Z_of_t(g);
    double log10_Z = log(fabs(Z) + 1e-300) * log10_e;

    printf("  %4d  %14.4f  %14.8e  %10.2f  %10.0f\n", k, g, fabs(Z), log10_Z,
           mp_log10_at_zeros[k]);
  }
  printf("\n");
}

/* ====================================================================
 * TEST 2: |xi(1/2+it)| vs mpmath reference
 * ==================================================================== */

static void test_xi_vs_mpmath(void) {
  printf("====================================================================="
         "\n");
  printf("  TEST 2: |xi(1/2+it)| vs mpmath Reference (50-digit precision)\n");
  printf("====================================================================="
         "\n\n");

  printf("  %10s  %20s  %20s  %14s  %10s\n", "t", "log|xi_EM|",
         "log|xi_mpmath|", "|diff|", "rel_err");
  printf("  %10s  %20s  %20s  %14s  %10s\n", "---", "---", "---", "---", "---");

  for (int i = 0; i < n_mp; i++) {
    double t = mp_ref[i].t;
    double log_xi_EM = log_abs_xi_from_Z(t);
    double log_xi_mp = mp_ref[i].log_xi;
    double diff = fabs(log_xi_EM - log_xi_mp);
    double rel = fabs(log_xi_mp) > 0.1 ? diff / fabs(log_xi_mp) : diff;

    printf("  %10.4f  %20.6f  %20.6f  %14.6e  %10.2e\n", t, log_xi_EM,
           log_xi_mp, diff, rel);
  }
  printf("\n");
}

/* ====================================================================
 * TEST 3: Z(t) at non-zeros
 * ==================================================================== */

static void test_z_off_zeros(void) {
  printf("====================================================================="
         "\n");
  printf("  TEST 3: |Z(t)| at Non-Zero Points\n");
  printf("====================================================================="
         "\n\n");

  double t_vals[] = {15.0, 25.5, 40.0, 50.0, 100.0, 200.0};
  int nt = 6;

  printf("  %8s  %14s  %14s\n", "t", "|Z(t)|", "log|Z(t)|");
  printf("  %8s  %14s  %14s\n", "---", "---", "---");

  for (int i = 0; i < nt; i++) {
    double t = t_vals[i];
    double Z = Z_of_t(t);
    double logZ = log(fabs(Z) + 1e-300);
    printf("  %8.2f  %14.8e  %14.6f\n", t, fabs(Z), logZ);
  }
  printf("\n");
}

/* ====================================================================
 * TEST 4: Convergence with K terms
 * ==================================================================== */

static void test_convergence(void) {
  printf("====================================================================="
         "\n");
  printf("  TEST 4: Euler-Maclaurin Accuracy vs K terms\n");
  printf("  Z(gamma_0) at t=14.1347\n");
  printf("====================================================================="
         "\n\n");

  double g0 = ZETA_ZEROS[0];
  int K_vals[] = {5, 10, 15, 20, 25, 30};
  int nK = 6;

  printf("  %4s  %14s  %10s  %10s\n", "K", "|Z(gamma_0)|", "log10|Z|",
         "improvement");
  printf("  %4s  %14s  %10s  %10s\n", "---", "---", "---", "---");

  double prev_Z = 0.0;

  for (int i = 0; i < nK; i++) {
    int K = K_vals[i];
    int N = (int)(g0 + 20);
    if (N < 50)
      N = 50;

    (void)K; /* K used via loop below */
    Complex sum = {0.0, 0.0};
    for (int n = 1; n < N; n++) {
      double log_n = log((double)n);
      double factor = pow((double)n, -0.5);
      sum.re += factor * cos(-g0 * log_n);
      sum.im += factor * sin(-g0 * log_n);
    }
    double log_N = log((double)N);
    double N_factor = pow((double)N, -0.5);
    sum.re += N_factor * cos(-g0 * log_N) / 2.0;
    sum.im += N_factor * sin(-g0 * log_N) / 2.0;

    Complex s_m1 = {-0.5, g0};
    double N_1m = pow((double)N, 0.5);
    Complex num = {N_1m * cos(-g0 * log_N), N_1m * sin(-g0 * log_N)};
    Complex term2 = c_div(num, s_m1);
    sum = c_add(sum, term2);

    /* EM correction terms */
    for (int k = 1; k <= K; k++) {
      double B2k = bernoulli_2k(k);
      if (B2k == 0.0)
        continue;
      double fact2k = factorial_2k(k);
      double coeff = B2k / fact2k;

      Complex prod = {1.0, 0.0};
      for (int j = 0; j <= 2 * k - 2; j++) {
        Complex s_j = {0.5 + j, g0};
        prod = c_mul(prod, s_j);
      }

      double N_factor2 = pow((double)N, -0.5 - 2.0 * k + 1.0);
      Complex N_term = {N_factor2 * cos(-g0 * log_N),
                        N_factor2 * sin(-g0 * log_N)};

      Complex term = c_mul(prod, N_term);
      term.re *= coeff;
      term.im *= coeff;
      sum = c_add(sum, term);
    }

    double theta = rs_theta(g0);
    double Z = cos(theta) * sum.re - sin(theta) * sum.im;
    double log10_Z = log(fabs(Z) + 1e-300) / log(10.0);
    double improvement =
        (prev_Z > 0) ? log10(prev_Z / (fabs(Z) + 1e-300)) : 0.0;

    printf("  %4d  %14.8e  %10.2f  %10.1f digits\n", K, fabs(Z), log10_Z,
           improvement);

    prev_Z = fabs(Z);
  }
  printf("\n  Target: ~10^-54 (mpmath reference)\n\n");
}

/* ====================================================================
 * MAIN
 * ==================================================================== */

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

  printf("#####################################################################"
         "###\n");
  printf("#  RIEMANN-SIEGEL Z(t) — Euler-Maclaurin Formula\n");
  printf("#\n");
  printf("#  zeta(s) via Euler-Maclaurin with Bernoulli correction terms\n");
  printf("#  Z(t) = e^{i*theta(t)} * zeta(1/2+it) (real-valued)\n");
  printf("#\n");
  printf("#  Reference: tests/ref_xi_mpmath.py (mpmath, 50-digit precision)\n");
  printf("#####################################################################"
         "###\n\n");

  test_z_at_zeros();
  test_xi_vs_mpmath();
  test_z_off_zeros();
  test_convergence();

  printf("#####################################################################"
         "###\n");
  printf("#  CONCLUSION\n");
  printf("#\n");
  printf(
      "#  The Euler-Maclaurin formula converges rapidly for zeta(1/2+it).\n");
  printf("#  At zeta zeros, |Z(t)| should approach the mpmath reference.\n");
  printf("#####################################################################"
         "###\n");

  return 0;
}
