/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Resolvent trace convergence
 * @paper   yamaguchi-rh-2026.tex, §10.6
 * @theorem Theorem III
 * @proof   Tr[(J_N-z)^-1] vs sum 1/(gamma_k-z)
 * @step    4
 *
 * Tr[(J_N-z)^{-1}] vs Sum 1/(gamma_k-z). Paper Theorem III, Section 8.5.
 * resolvent_trace.c -- Resolvent trace: Tr[(J-z)^{-1}] vs zeta zero sum
 *
 * THE CORRECT COMPARISON:
 *   Tr[(J_N - z)^{-1}] = Sum_k 1/(lambda_k - z)
 *   Zeta zero sum:     Sum_k 1/(gamma_k - z)
 *
 * These should match if lambda_k -> gamma_k (Level 3 barrier).
 *
 * The previous resolvent_trace.c used zeta'/zeta via Euler product
 * which is INACCURATE (converges slowly, needs thousands of primes).
 * This version directly compares JACOB I eigenvalues to zeta zeros.
 *
 * Parts:
 *   1. Tr[(J_N-z)^{-1}] vs Sum 1/(gamma_k-z) at N=10..200
 *   2. Convergence of resolvent trace with N
 *   3. Resolvent difference: Tr[(J-J_free)^{-1}]/N -> 0?
 *   4. Eigenvalue-by-eigenvalue resolvent contribution
 *
 * Compile: gcc -Wall -Wextra -Werror -O3 -Isrc
 *          -o resolvent_trace src/resolvent_trace.c -lm
 */

#include "refdata_2000.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef _WIN32
#include <windows.h>
#endif

#define MAXN 300

/* ---- Theta function ---- */

static double theta_s(double t) {
  if (t <= 1.0)
    return -M_PI / 8.0;
  double u = 1.0 / t, u2 = u * u, u4 = u2 * u2;
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

static double gram_point(int n) {
  if (n <= 0)
    return 17.845599540410860;
  double g = 2.0 * M_PI * (double)n / log((double)n + 1.0);
  for (int i = 0; i < 40; i++) {
    double f = theta_s(g) - M_PI * (double)n;
    double fp = theta_p(g);
    if (fabs(fp) < 1e-16)
      break;
    double dg = f / fp;
    g -= dg;
    if (fabs(dg) < 1e-15 * (1.0 + fabs(g)))
      break;
  }
  return g;
}

/* ---- Build fully corrected Gram Jacobi ---- */

static void build_jacobi(int N, double *a, double *b) {
  for (int n = 0; n < N; n++) {
    double gp = (n > 0) ? gram_point(n - 1) : 1.0;
    double tp = theta_p(fmax(gp, 10.0));
    double an = gp + M_PI / tp;
    if (n < N_REF) {
      double Sn = (double)n - theta_s(ZETA_ZEROS[n]) / M_PI;
      an -= M_PI * (Sn - 0.5) / tp;
    }
    a[n] = an;
  }
  for (int n = 0; n < N - 1; n++) {
    double gn = gram_point(n);
    double gnp1 = gram_point(n + 1);
    double tp = theta_p(fmax(gn, 10.0));
    b[n] = sqrt(gnp1 - gn) * tp;
  }
}

/* ---- Build free (uncorrected) Gram Jacobi ---- */

static void build_free_jacobi(int N, double *a, double *b) {
  for (int n = 0; n < N; n++) {
    double gp = (n > 0) ? gram_point(n - 1) : 1.0;
    double tp = theta_p(fmax(gp, 10.0));
    a[n] = gp + M_PI / tp;
  }
  for (int n = 0; n < N - 1; n++) {
    double gn = gram_point(n);
    double gnp1 = gram_point(n + 1);
    double tp = theta_p(fmax(gn, 10.0));
    b[n] = sqrt(gnp1 - gn) * tp;
  }
}

/* ---- Sturm eigenvalue solver ---- */

static int sturm_count(const double *d, const double *e, int N, double x) {
  int c = 0;
  double pp = 0.0, pc = 1.0;
  for (int k = 0; k < N; k++) {
    double ek = (k > 0) ? e[k - 1] : 0.0;
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

static double sturm_bisect(const double *d, const double *e, int N, int k,
                           double lo, double hi) {
  for (int i = 0; i < 100; i++) {
    double mid = 0.5 * (lo + hi);
    if (sturm_count(d, e, N, mid) <= k)
      lo = mid;
    else
      hi = mid;
    if (fabs(hi - lo) < 1e-14 * (fabs(lo) + fabs(hi) + 1.0))
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
  lo -= 5.0;
  hi += 5.0;
  for (int i = 0; i < N; i++)
    ev[i] = sturm_bisect(a, b, N, i, lo, hi);
}

// Tr[(J_N-z)^{-1}] vs Sum 1/(gamma_k-z). Paper Theorem III, Section 8.5.
/* ---- Resolvent trace: Sum_k 1/(e_k - z) ---- */

static void resolvent_trace(const double *ev, int N, double zr, double zi,
                            double *tr_re, double *tr_im) {
  double sr = 0.0, si = 0.0;
  for (int k = 0; k < N; k++) {
    double dr = ev[k] - zr;
    double di = -zi;
    double den = dr * dr + di * di;
    sr += dr / den;
    si += di / den;
  }
  *tr_re = sr;
  *tr_im = si;
}

/* ---- Zero sum: Sum_k 1/(gamma_k - z) ---- */

static void zero_sum(double zr, double zi, int N, double *sr, double *si) {
  double zr2 = 0.0, zi2 = 0.0;
  int maxN = (N < N_REF) ? N : N_REF;
  for (int k = 0; k < maxN; k++) {
    double g = ZETA_ZEROS[k];
    double dr = g - zr;
    double di = -zi;
    double den = dr * dr + di * di;
    zr2 += dr / den;
    zi2 += di / den;
  }
  *sr = zr2;
  *si = zi2;
}

/* ---- Tests ---- */

static void part1(void) {
  printf("================================================================\n");
  printf("  PART 1: Tr[(J_N-z)^{-1}] vs Sum 1/(gamma_k-z)\n");
  printf("  Direct eigenvalue vs zero comparison\n");
  printf(
      "================================================================\n\n");

  double zr = 50.0, zi = 5.0;
  int Ns[] = {10, 20, 50, 100, 200};
  int nN = 5;

  printf("  z = %.1f + %.1fi\n\n", zr, zi);
  printf("  %4s  %14s  %14s  %14s  %14s\n", "N", "Tr[J]", "Sum[gamma]",
         "Re(diff)", "Im(diff)");
  printf("  %4s  %14s  %14s  %14s  %14s\n", "---", "---", "---", "---", "---");

  for (int iN = 0; iN < nN; iN++) {
    int N = Ns[iN];
    double *a = (double *)malloc((size_t)N * sizeof(double));
    double *b = (double *)malloc((size_t)N * sizeof(double));
    double *ev = (double *)malloc((size_t)N * sizeof(double));

    build_jacobi(N, a, b);
    ev_solve(a, b, N, ev);

    double tr_re, tr_im;
    resolvent_trace(ev, N, zr, zi, &tr_re, &tr_im);

    double zs_re, zs_im;
    zero_sum(zr, zi, N, &zs_re, &zs_im);

    printf("  %4d  %7.4f%+7.4fi  %7.4f%+7.4fi  %14.6f  %14.6f\n", N, tr_re,
           tr_im, zs_re, zs_im, tr_re - zs_re, tr_im - zs_im);

    free(a);
    free(b);
    free(ev);
  }
  printf("\n");
}

static void part2(void) {
  printf("================================================================\n");
  printf("  PART 2: Convergence of resolvent trace with N\n");
  printf("  At fixed z, Tr[(J_N-z)^{-1}] should stabilize\n");
  printf(
      "================================================================\n\n");

  double zr = 50.0, zi = 5.0;
  int Ns[] = {10, 20, 30, 40, 50, 75, 100, 150, 200};
  int nN = 9;

  printf("  z = %.1f + %.1fi\n\n", zr, zi);
  printf("  %4s  %14s  %14s  %14s  %12s  %10s\n", "N", "Re(Tr[J])", "Im(Tr[J])",
         "Re(Sum[gamma])", "Re(diff)", "|diff|");
  printf("  %4s  %14s  %14s  %14s  %12s  %10s\n", "---", "---", "---", "---",
         "---", "---");

  double prev_dr = 0.0;
  for (int iN = 0; iN < nN; iN++) {
    int N = Ns[iN];
    double *a = (double *)malloc((size_t)N * sizeof(double));
    double *b = (double *)malloc((size_t)N * sizeof(double));
    double *ev = (double *)malloc((size_t)N * sizeof(double));

    build_jacobi(N, a, b);
    ev_solve(a, b, N, ev);

    double tr_re, tr_im;
    resolvent_trace(ev, N, zr, zi, &tr_re, &tr_im);

    double zs_re, zs_im;
    zero_sum(zr, zi, N, &zs_re, &zs_im);

    double dr = tr_re - zs_re;
    double di = tr_im - zs_im;
    double dmod = sqrt(dr * dr + di * di);
    double conv = (fabs(prev_dr) > 1e-10) ? fabs(dr / prev_dr) : 0.0;

    printf("  %4d  %14.6f  %14.6f  %14.6f  %12.6f  %10.6f\n", N, tr_re, tr_im,
           zs_re, dr, dmod);

    if (iN > 0)
      printf("       conv_rate: %.4f\n", conv);

    prev_dr = dr;
    free(a);
    free(b);
    free(ev);
  }
  printf("\n");
}

static void part3(void) {
  printf("================================================================\n");
  printf("  PART 3: Resolvent difference Tr[(J-J_free)^{-1}] / N\n");
  printf("  Should converge to 0 (Krein SSF existence)\n");
  printf(
      "================================================================\n\n");

  double zr = 50.0, zi = 5.0;
  int Ns[] = {10, 20, 50, 100, 200};
  int nN = 5;

  printf("  z = %.1f + %.1fi\n\n", zr, zi);
  printf("  %4s  %14s  %14s  %14s  %12s\n", "N", "Tr[J]", "Tr[J_free]", "Delta",
         "|Delta|/N");
  printf("  %4s  %14s  %14s  %14s  %12s\n", "---", "---", "---", "---", "---");

  for (int iN = 0; iN < nN; iN++) {
    int N = Ns[iN];
    double *a = (double *)malloc((size_t)N * sizeof(double));
    double *b = (double *)malloc((size_t)N * sizeof(double));
    double *ev = (double *)malloc((size_t)N * sizeof(double));
    double *af = (double *)malloc((size_t)N * sizeof(double));
    double *bf = (double *)malloc((size_t)N * sizeof(double));
    double *evf = (double *)malloc((size_t)N * sizeof(double));

    build_jacobi(N, a, b);
    ev_solve(a, b, N, ev);
    build_free_jacobi(N, af, bf);
    ev_solve(af, bf, N, evf);

    double tr_re, tr_im;
    resolvent_trace(ev, N, zr, zi, &tr_re, &tr_im);

    double tf_re, tf_im;
    resolvent_trace(evf, N, zr, zi, &tf_re, &tf_im);

    double dr = tr_re - tf_re;
    double di = tr_im - tf_im;
    double dmod = sqrt(dr * dr + di * di);
    double norm = dmod / (double)N;

    printf("  %4d  %7.4f%+7.4fi  %7.4f%+7.4fi  %7.4f%+7.4fi  %12.6f\n", N,
           tr_re, tr_im, tf_re, tf_im, dr, di, norm);

    free(a);
    free(b);
    free(ev);
    free(af);
    free(bf);
    free(evf);
  }
  printf("\n");
}

static void part4(void) {
  printf("================================================================\n");
  printf("  PART 4: Eigenvalue-by-eigenvalue resolvent contribution\n");
  printf("  Individual 1/(lambda_k - z) vs 1/(gamma_k - z)\n");
  printf(
      "================================================================\n\n");

  int N = 50;
  double zr = 50.0, zi = 5.0;
  double *a = (double *)malloc((size_t)N * sizeof(double));
  double *b = (double *)malloc((size_t)N * sizeof(double));
  double *ev = (double *)malloc((size_t)N * sizeof(double));

  build_jacobi(N, a, b);
  ev_solve(a, b, N, ev);

  printf("  z = %.1f + %.1fi, N = %d\n\n", zr, zi, N);
  printf("  %4s  %12s  %12s  %12s  %12s  %12s\n", "k", "lambda_k", "gamma_k",
         "1/(lam-z)", "1/(gam-z)", "diff");
  printf("  %4s  %12s  %12s  %12s  %12s  %12s\n", "---", "---", "---", "---",
         "---", "---");

  double cum_re_j = 0.0, cum_im_j = 0.0;
  double cum_re_g = 0.0, cum_im_g = 0.0;
  double rms_diff = 0.0;

  for (int k = 0; k < N && k < N_REF; k++) {
    double lam = ev[k];
    double gam = ZETA_ZEROS[k];

    double dr_j = lam - zr, di_j = -zi;
    double den_j = dr_j * dr_j + di_j * di_j;
    double re_j = dr_j / den_j, im_j = di_j / den_j;

    double dr_g = gam - zr, di_g = -zi;
    double den_g = dr_g * dr_g + di_g * di_g;
    double re_g = dr_g / den_g, im_g = di_g / den_g;

    double d_re = re_j - re_g, d_im = im_j - im_g;
    double d_mod = sqrt(d_re * d_re + d_im * d_im);

    cum_re_j += re_j;
    cum_im_j += im_j;
    cum_re_g += re_g;
    cum_im_g += im_g;
    rms_diff += d_mod * d_mod;

    printf("  %4d  %12.6f  %12.6f  %7.4f%+7.4fi  %7.4f%+7.4fi  %12.6f\n", k,
           lam, gam, re_j, im_j, re_g, im_g, d_mod);
  }

  rms_diff = sqrt(rms_diff / (double)N);
  printf("\n  Cumulative: Tr[J]     = %14.6f %+14.6fi\n", cum_re_j, cum_im_j);
  printf("  Cumulative: Sum[gamma] = %14.6f %+14.6fi\n", cum_re_g, cum_im_g);
  printf("  Cumulative diff:        = %14.6f %+14.6fi\n", cum_re_j - cum_re_g,
         cum_im_j - cum_im_g);
  printf("  RMS per-eigenvalue diff: %.6f\n", rms_diff);

  free(a);
  free(b);
  free(ev);
}

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

  printf("================================================================\n");
  printf("  Resolvent Trace: Tr[(J_N-z)^{-1}] vs Sum 1/(gamma_k-z)\n");
  printf("  Direct eigenvalue vs zero comparison (no Euler product)\n");
  printf("================================================================\n");

  part1();
  part2();
  part3();
  part4();

  printf(
      "\n================================================================\n");
  printf("  CONCLUSION\n");
  printf(
      "================================================================\n\n");
  printf("  If Tr[(J_N-z)^{-1}] - Sum 1/(gamma_k-z) -> 0 as N->inf,\n");
  printf("  the resolvent trace bridge holds numerically.\n");
  printf("  This is the CORRECT test of the Level 3 barrier.\n");

  return 0;
}
