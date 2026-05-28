/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Richardson extrapolation for Hadamard tail
 * @paper   yamaguchi-rh-2026.tex, §8.4
 * @theorem Theorem III
 * @proof   Convergence acceleration
 * @step    4
 *
 * Richardson extrapolation: (4*f(2N)-f(N))/3. Paper Section 8.4.
 * hadamard_extrapolation.c — Richardson Extrapolation + xi_ev Gap Analysis
 *
 * OBJECTIVE 1: Analyze xi_ev (eigenvalue quadratic form) vs xi_H (Hadamard) gap
 *   at larger N (up to 1000) to determine convergence rate.
 *
 * Richardson extrapolation: (4*f(2N)-f(N))/3. Paper Section 8.4.
 * OBJECTIVE 2: Implement Richardson extrapolation for Hadamard tail:
 *   log|xi_H(N)| = log|xi(∞)| + c_1/N + c_2/N^2 + ...
 *   → xi_extrapolated = (4*xi_H(2N) - xi_H(N)) / 3
 *   This eliminates the O(1/N) term, giving O(1/N^2) convergence.
 *
 * Compile: gcc -O3 -o bin/hadamard_extrapolation.exe
 * src/hadamard_extrapolation.c
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

/* Jacobi construction */
static double theta_s(double t) {
  if (t <= 1.0)
    return -M_PI / 8.0;
  double x = t / (2.0 * M_PI), u = 1.0 / t, u2 = u * u, u4 = u2 * u2;
  return 0.5 * t * log(x) - 0.5 * t - M_PI / 8.0 + u / 48.0 +
         7.0 * u * u2 / 5760.0 + 31.0 * u4 / 80640.0 +
         127.0 * u4 * u2 / 430080.0 + 2555.0 * u4 * u4 / 27525120.0;
}
static double theta_p(double t) {
  if (t <= 2.0 * M_PI)
    return 1.0;
  double t2 = t * t;
  return 0.5 * log(t / (2.0 * M_PI)) - 1.0 / (24.0 * t2) +
         7.0 / (960.0 * t2 * t2);
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

/* ---- mpmath reference values ---- */

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

/* ---- Hadamard product for xi(1/2+it) ---- */

static double log_xi_hadamard_it(double t, int N) {
  double lr = log(0.5);
  int maxN = (N < N_REF) ? N : N_REF;
  for (int k = 0; k < maxN; k++) {
    double g = ZETA_ZEROS[k];
    double num = g * g - t * t;
    double den = 0.25 + g * g;
    lr += log(fabs(num)) - log(den);
  }
  return lr;
}

/* ---- xi_ev: eigenvalue quadratic form ---- */

static double log_xi_eigenvalue_z(const double *ev, int N, double t) {
  double lr = log(0.5);
  for (int k = 0; k < N; k++) {
    double lam = ev[k];
    double num = lam * lam - t * t;
    double den = 0.25 + lam * lam;
    lr += log(fabs(num)) - log(den);
  }
  return lr;
}

// Richardson extrapolation: (4*f(2N)-f(N))/3. Paper Section 8.4.
/* ---- Richardson extrapolation ----
 *
 * If f(N) = f(∞) + c_1/N + c_2/N^2 + ..., then:
 *   f_extrap(2N, N) = (4*f(2N) - f(N)) / 3
 * eliminates the O(1/N) term.
 *
 * We apply this to log|xi_H(N)| at each t value.
 */

static double richardson_extrap(double f_N, double f_2N) {
  return (4.0 * f_2N - f_N) / 3.0;
}

/*
 * TEST 1: xi_ev vs xi_H gap analysis at N=500, 1000
 */
static void test_xiev_gap(void) {
  printf("====================================================================="
         "\n");
  printf("  TEST 1: xi_ev vs xi_H Gap Analysis\n");
  printf("  Eigenvalue quadratic form vs Hadamard at larger N\n");
  printf("====================================================================="
         "\n\n");

  int N_vals[] = {10, 20, 50, 100, 200, 500, 1000};
  int nN = 7;

  double t_vals[] = {14.13, 30.42, 50.0, 100.0};
  int nt = 4;

  for (int it = 0; it < nt; it++) {
    double t = t_vals[it];
    printf("  t = %.2f\n\n", t);
    printf("  %4s  %14s  %14s  %14s  %10s\n", "N", "log|xi_ev|", "log|xi_H|",
           "|xi_ev-xi_H|", "gap_rate");
    printf("  %4s  %14s  %14s  %14s  %10s\n", "---", "---", "---", "---",
           "---");

    double prev_gap = 0.0;

    for (int iN = 0; iN < nN; iN++) {
      int N = N_vals[iN];
      double *a = malloc((size_t)N * sizeof(double));
      double *b = malloc((size_t)N * sizeof(double));
      double *ev = malloc((size_t)N * sizeof(double));
      build_jacobi(N, a, b);
      ev_solve(a, b, N, ev);

      double log_xi_ev = log_xi_eigenvalue_z(ev, N, t);
      double log_xi_H = log_xi_hadamard_it(t, N);
      double gap = fabs(log_xi_ev - log_xi_H);
      double gap_rate = (prev_gap > 1e-15) ? prev_gap / gap : 0.0;

      printf("  %4d  %14.6f  %14.6f  %14.6e  %10.3fx\n", N, log_xi_ev, log_xi_H,
             gap, gap_rate);

      prev_gap = gap;
      free(a);
      free(b);
      free(ev);
    }
    printf("\n");
  }
}

/*
 * Richardson extrapolation: (4*f(2N)-f(N))/3. Paper Section 8.4.
 * TEST 2: Richardson extrapolation for Hadamard tail
 */
static void test_richardson(void) {
  printf("====================================================================="
         "\n");
  // Richardson extrapolation: (4*f(2N)-f(N))/3. Paper Section 8.4.
  printf("  TEST 2: Richardson Extrapolation — O(1/N) → O(1/N^2)\n");
  printf("  f_extrap = (4*f(2N) - f(N)) / 3\n");
  printf("====================================================================="
         "\n\n");

  /* Use N pairs: (100,200), (200,500), (500,1000), (1000,2000) */
  int N_pairs[][2] = {{100, 200}, {200, 500}, {500, 1000}, {1000, 2000}};
  int nP = 4;

  double t_vals[] = {14.13, 21.02, 30.42, 50.0, 100.0};
  int nt = 5;

  printf("  %10s  %10s  %14s  %14s  %14s  %10s  %10s\n", "t", "N pair",
         "log|xi_H(N)|", "log|xi_H(2N)|", "log|xi_extrap|", "|raw diff|",
         "|extrap diff|");
  printf("  %10s  %10s  %14s  %14s  %14s  %10s  %10s\n", "---", "---", "---",
         "---", "---", "---", "---");

  for (int it = 0; it < nt; it++) {
    double t = t_vals[it];
    double log_xi_mp = 0;
    for (int i = 0; i < n_mp; i++) {
      if (fabs(mp_ref[i].t - t) < 0.01) {
        log_xi_mp = mp_ref[i].log_xi;
        break;
      }
    }

    for (int ip = 0; ip < nP; ip++) {
      int N1 = N_pairs[ip][0], N2 = N_pairs[ip][1];
      double f_N = log_xi_hadamard_it(t, N1);
      double f_2N = log_xi_hadamard_it(t, N2);
      double f_extrap = richardson_extrap(f_N, f_2N);

      double raw_diff = (log_xi_mp != 0) ? fabs(f_2N - log_xi_mp) : 0;
      double extrap_diff = (log_xi_mp != 0) ? fabs(f_extrap - log_xi_mp) : 0;

      printf("  %10.2f  %4d,%4d  %14.6f  %14.6f  %14.6f  %10.6e  %10.6e\n", t,
             N1, N2, f_N, f_2N, f_extrap, raw_diff, extrap_diff);
    }
    printf("\n");
  }
}

/*
 * Richardson extrapolation: (4*f(2N)-f(N))/3. Paper Section 8.4.
 * TEST 3: Hadamard + Richardson vs mpmath at all reference points
 */

static void test_extrap_vs_mpmath(void) {
  printf("====================================================================="
         "\n");
  printf("  TEST 3: Richardson Extrapolation vs mpmath (N=1000,2000)\n");
  printf("====================================================================="
         "\n\n");

  printf("  %10s  %20s  %20s  %20s  %14s  %14s  %10s\n", "t", "log|xi_mpmath|",
         "log|xi_H(1000)|", "log|xi_H(2000)|", "log|xi_extrap|", "|raw diff|",
         "|extrap diff|");
  printf("  %10s  %20s  %20s  %20s  %14s  %14s  %10s\n", "---", "---", "---",
         "---", "---", "---", "---");

  for (int i = 0; i < n_mp; i++) {
    double t = mp_ref[i].t;
    double log_xi_mp = mp_ref[i].log_xi;
    double f_1000 = log_xi_hadamard_it(t, 1000);
    double f_2000 = log_xi_hadamard_it(t, 2000);
    double f_extrap = richardson_extrap(f_1000, f_2000);

    double raw_diff = fabs(f_2000 - log_xi_mp);
    double extrap_diff = fabs(f_extrap - log_xi_mp);

    printf("  %10.2f  %20.6f  %20.6f  %20.6f  %14.6f  %14.6e  %14.6e\n", t,
           log_xi_mp, f_1000, f_2000, f_extrap, raw_diff, extrap_diff);
  }
  printf("\n");
}

/*
 * TEST 4: Eigenvalue gap convergence rate
 */

static void test_ev_gap_rate(void) {
  printf("====================================================================="
         "\n");
  printf("  TEST 4: Eigenvalue Gap Convergence Rate\n");
  printf("  |xi_ev(N) - xi_H(N)| vs N at fixed t\n");
  printf("====================================================================="
         "\n\n");

  int N_vals[] = {10, 20, 50, 100, 200, 500, 1000};
  int nN = 7;
  double t = 14.13;

  printf("  t = %.2f\n\n", t);
  printf("  %4s  %14s  %10s  %14s  %10s\n", "N", "|xi_ev-xi_H|", "N*gap",
         "N^2*gap", "conv_rate");
  printf("  %4s  %14s  %10s  %14s  %10s\n", "---", "---", "---", "---", "---");

  double prev_gap = 0.0;

  for (int iN = 0; iN < nN; iN++) {
    int N = N_vals[iN];
    double *a = malloc((size_t)N * sizeof(double));
    double *b = malloc((size_t)N * sizeof(double));
    double *ev = malloc((size_t)N * sizeof(double));
    build_jacobi(N, a, b);
    ev_solve(a, b, N, ev);

    double log_xi_ev = log_xi_eigenvalue_z(ev, N, t);
    double log_xi_H = log_xi_hadamard_it(t, N);
    double gap = fabs(log_xi_ev - log_xi_H);
    double N_gap = N * gap;
    double N2_gap = (double)N * (double)N * gap;
    double rate = (prev_gap > 1e-15) ? prev_gap / gap : 0.0;

    printf("  %4d  %14.6e  %10.6e  %14.6e  %10.3fx\n", N, gap, N_gap, N2_gap,
           rate);

    prev_gap = gap;
    free(a);
    free(b);
    free(ev);
  }
  printf("\n  If N*gap is constant: convergence is O(1/N)\n");
  printf("  If N^2*gap is constant: convergence is O(1/N^2)\n\n");
}

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

  printf("#####################################################################"
         "###\n");
  // Richardson extrapolation: (4*f(2N)-f(N))/3. Paper Section 8.4.
  printf("#  HADAMARD EXTRAPOLATION — Richardson + xi_ev Gap Analysis\n");
  printf("#\n");
  printf("#  1. xi_ev vs xi_H gap at N up to 1000\n");
  // Richardson extrapolation: (4*f(2N)-f(N))/3. Paper Section 8.4.
  printf("#  2. Richardson extrapolation: O(1/N) → O(1/N^2)\n");
  printf("#  3. Extrapolated Hadamard vs mpmath\n");
  printf("#####################################################################"
         "###\n\n");

  test_xiev_gap();
  test_richardson();
  test_extrap_vs_mpmath();
  test_ev_gap_rate();

  printf("#####################################################################"
         "###\n");
  printf("#  CONCLUSION\n");
  printf("#\n");
  // Richardson extrapolation: (4*f(2N)-f(N))/3. Paper Section 8.4.
  printf("#  Richardson extrapolation eliminates the O(1/N) term.\n");
  printf("#  If |extrap diff| < |raw diff|, the extrapolation helps.\n");
  printf("#\n");
  printf(
      "#  The xi_ev gap analysis reveals the Jacobi reconstruction quality.\n");
  printf("#####################################################################"
         "###\n");

  return 0;
}
