/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   m-function bridge: m(z) vs zeta'/zeta
 * @paper   yamaguchi-rh-2026.tex, §10.6
 * @theorem Theorem III
 * @proof   Continued fraction + Krein SSF + entry DFT
 * @step    4
 *
 * mfunction_bridge.c — Three critical computations for the RH proof chain
 *
 * Part 1: m(z) vs zeta'/zeta direct comparison
 * Part 2: Krein SSF pointwise test
 * Part 3: Entry analytic formula — prime coefficient extraction via DFT
 *
 * Compile: gcc -Wall -Wextra -Wconversion -Wshadow -Werror -O3
 *          -fno-strict-aliasing -fno-peel-loops -fno-unswitch-loops
 *          -I. -o mfunction_bridge.exe src/mfunction_bridge.c -lm
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

#define MAXN 50
#define NPRIMES200 200

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

static const int PRIMES200[NPRIMES200] = {
    2,    3,    5,    7,    11,   13,   17,   19,   23,   29,   31,   37,
    41,   43,   47,   53,   59,   61,   67,   71,   73,   79,   83,   89,
    97,   101,  103,  107,  109,  113,  127,  131,  137,  139,  149,  151,
    157,  163,  167,  173,  179,  181,  191,  193,  197,  199,  211,  223,
    227,  229,  233,  239,  241,  251,  257,  263,  269,  271,  277,  281,
    283,  293,  307,  311,  313,  317,  331,  337,  347,  349,  353,  359,
    367,  373,  379,  383,  389,  397,  401,  409,  419,  421,  431,  433,
    439,  443,  449,  457,  461,  463,  467,  479,  487,  491,  499,  503,
    509,  521,  523,  541,  547,  557,  563,  569,  571,  577,  587,  593,
    599,  601,  607,  613,  617,  619,  631,  641,  643,  647,  653,  659,
    661,  673,  677,  683,  691,  701,  709,  719,  727,  733,  739,  743,
    751,  757,  761,  769,  773,  787,  797,  809,  811,  821,  823,  827,
    829,  839,  853,  857,  859,  863,  877,  881,  883,  887,  907,  911,
    919,  929,  937,  941,  947,  953,  967,  971,  977,  983,  991,  997,
    1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
    1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163,
    1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223};

static int deboor(const double *lam, const double *mu, int N, double *a,
                  double *b) {
  for (int k = 0; k < N - 1; k++)
    if (!(lam[k] < mu[k] && mu[k] < lam[k + 1]))
      return -1;
  double w[MAXN];
  double ws = 0.0;
  for (int k = 0; k < N; k++) {
    double nm = 1.0;
    for (int j = 0; j < N - 1; j++)
      nm *= lam[k] - mu[j];
    double dn = 1.0;
    for (int j = 0; j < N; j++)
      if (j != k)
        dn *= lam[k] - lam[j];
    w[k] = nm / dn;
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

static double theta_s(double t) {
  if (t <= 1.0)
    return -M_PI / 8.0;
  double x = t / (2.0 * M_PI);
  double u = 1.0 / t;
  return 0.5 * t * log(x) - 0.5 * t - M_PI / 8.0 + u / 48.0 +
         7.0 * u * u * u / 5760.0;
}

static double theta_p(double t) {
  if (t <= 2.0 * M_PI)
    return 1.0;
  return 0.5 * log(t / (2.0 * M_PI)) - 1.0 / (24.0 * t * t);
}

__attribute__((unused)) static double gram_pt(int n) {
  double g = (n < 2) ? 17.8456 : 2.0 * M_PI * (double)n / log((double)n + 1.0);
  for (int i = 0; i < 80; i++) {
    double f = theta_s(g) - M_PI * (double)n;
    double d = f / theta_p(g);
    if (fabs(d) < 1e-14 * (1.0 + fabs(g)))
      break;
    g -= d;
  }
  return g;
}

static Cpx zeta_prime_over_zeta(double E, int np) {
  Cpx sum = cpx_make(0.0, 0.0);
  for (int pi = 0; pi < np; pi++) {
    double p = (double)PRIMES200[pi];
    double lp = log(p);
    Cpx phase = cpx_p_ik(p, E);
    Cpx ps = cpx_scale(sqrt(p), phase);
    Cpx denom = cpx_sub(ps, cpx_make(1.0, 0.0));
    if (cpx_abs(denom) < 1e-10)
      continue;
    Cpx term = cpx_inv(denom);
    term = cpx_scale(-lp, term);
    sum = cpx_add(sum, term);
  }
  return sum;
}

static Cpx m_function(const double *a_arr, const double *b_arr, int N,
                      double zr, double zi) {
  double cr = a_arr[N - 1] - zr;
  double ci = -zi;
  double dn = cr * cr + ci * ci;
  if (dn < 1e-30)
    return cpx_make(1e15, 0.0);
  double mr = cr / dn;
  double mi = -ci / dn;
  for (int k = N - 2; k >= 0; k--) {
    double b2 = b_arr[k] * b_arr[k];
    double vr = a_arr[k] - zr - b2 * mr;
    double vi = -zi - b2 * mi;
    dn = vr * vr + vi * vi;
    if (dn < 1e-30)
      return cpx_make(1e15, 0.0);
    mr = vr / dn;
    mi = -vi / dn;
  }
  return cpx_make(mr, mi);
}

static int sturm_count(const double *d, const double *e, int N, double x) {
  int c = 0;
  double pp = 0.0, pc = 1.0;
  for (int k = 0; k < N; k++) {
    double ek = (k > 0) ? e[k - 1] : 0.0;
    double pn = (d[k] - x) * pc - ek * ek * pp;
    if (fabs(pn) > 1e100)
      pn = (pn > 0.0) ? 1e100 : -1e100;
    if (pc * pn < 0.0)
      c++;
    pp = pc;
    pc = pn;
  }
  return c;
}

static double arg_zeta_half(double E, int np) {
  double total_arg = 0.0;
  for (int pi = 0; pi < np; pi++) {
    double p = (double)PRIMES200[pi];
    Cpx phase = cpx_p_ik(p, E);
    double ps = sqrt(p);
    double re = 1.0 - phase.re / ps;
    double im = phase.im / ps;
    total_arg += atan2(im, re);
  }
  return -total_arg;
}

static double resolvent_trace_imag(const double *a_arr, const double *b_arr,
                                   int N, double E, double eps) {
  double tr_im = 0.0;
  for (int k = 0; k < N; k++) {
    double dr = a_arr[k] - E;
    double di = -eps;
    double d = dr * dr + di * di;
    tr_im += di / d;
  }
  for (int k = 0; k < N - 1; k++) {
    double b2 = b_arr[k] * b_arr[k];
    double re_sum = 0.0;
    double im_sum = 0.0;
    Cpx m_cur = cpx_inv(cpx_make(a_arr[N - 1] - E, -eps));
    for (int j = N - 2; j >= 0; j--) {
      re_sum += m_cur.re;
      im_sum += m_cur.im;
      if (j == k)
        break;
      double vr = a_arr[j] - E - b_arr[j] * b_arr[j] * m_cur.re;
      double vi = -eps - b_arr[j] * b_arr[j] * m_cur.im;
      double dn = vr * vr + vi * vi;
      if (dn < 1e-30)
        break;
      m_cur = cpx_make(vr / dn, -vi / dn);
    }
    tr_im += 2.0 * b2 * im_sum;
    (void)re_sum;
  }
  return tr_im;
}

static void part1(void) {
  printf("\n");
  printf("================================================================\n");
  printf("  PART 1: m(z) vs zeta'/zeta direct comparison\n");
  printf(
      "================================================================\n\n");

  int N = 30;
  double lam[MAXN], mu_arr[MAXN];
  for (int k = 0; k < N; k++)
    lam[k] = ZZ[k];
  for (int k = 0; k < N - 1; k++)
    mu_arr[k] = 0.5 * (ZZ[k] + ZZ[k + 1]);

  double a_arr[MAXN], b_arr[MAXN];
  if (deboor(lam, mu_arr, N, a_arr, b_arr) != 0) {
    printf("  dBG failed for N=%d\n", N);
    return;
  }

  int np = NPRIMES200;

  printf("  m(z) = 1/(a_0 - z - b_0^2/(a_1 - z - b_1^2/(.../(a_{N-1}-z))))\n");
  printf("  zeta'/zeta(1/2+iE) = -Sum_p ln(p)/(p^s - 1), p^s = sqrt(p)*e^{iE "
         "ln p}\n");
  printf("  Using %d primes, N=%d dBG entries\n\n", np, N);

  printf("  1a: At zeta zeros gamma_k (z = gamma_k + i*eps)\n");
  printf("  Both should blow up (residues at poles)\n\n");
  printf("  %4s  %10s  %5s  %12s  %12s  %12s  %12s  %12s\n", "k", "E", "eps",
         "Re(m)", "Im(m)", "|m|", "|z'/z|", "|Delta|");
  printf("  %4s  %10s  %5s  %12s  %12s  %12s  %12s  %12s\n", "---", "---",
         "---", "---", "---", "---", "---", "---");

  double eps_vals[] = {0.1, 0.01, 0.001};
  int n_eps = 3;

  for (int k = 0; k < 30; k++) {
    double E = ZZ[k];
    for (int ei = 0; ei < n_eps; ei++) {
      double eps = eps_vals[ei];
      Cpx mz = m_function(a_arr, b_arr, N, E, eps);
      Cpx zpz = zeta_prime_over_zeta(E, np);
      double delta = cpx_abs(cpx_sub(mz, zpz));
      printf("  %4d  %10.4f  %5.3f  %12.4f  %12.4f  %12.4f  %12.4f  %12.4f\n",
             k, E, eps, mz.re, mz.im, cpx_abs(mz), cpx_abs(zpz), delta);
    }
  }

  printf("\n  1b: At midpoints between zeros (regular points)\n\n");
  printf("  %4s  %10s  %5s  %12s  %12s  %12s  %12s  %12s\n", "k", "E_mid",
         "eps", "Re(m)", "Im(m)", "|m|", "|z'/z|", "|Delta|");
  printf("  %4s  %10s  %5s  %12s  %12s  %12s  %12s  %12s\n", "---", "---",
         "---", "---", "---", "---", "---", "---");

  for (int k = 0; k < 29; k++) {
    double E_mid = 0.5 * (ZZ[k] + ZZ[k + 1]);
    double eps = 0.1;
    Cpx mz = m_function(a_arr, b_arr, N, E_mid, eps);
    Cpx zpz = zeta_prime_over_zeta(E_mid, np);
    double delta = cpx_abs(cpx_sub(mz, zpz));
    printf("  %4d  %10.4f  %5.3f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f\n", k,
           E_mid, eps, mz.re, mz.im, cpx_abs(mz), cpx_abs(zpz), delta);
  }

  printf("\n  1c: Convergence with N at fixed z=50+0.1i\n\n");
  printf("  %6s  %12s  %12s  %12s  %12s  %12s\n", "N", "Re(m)", "Im(m)", "|m|",
         "|z'/z|", "|Delta|");
  printf("  %6s  %12s  %12s  %12s  %12s  %12s\n", "---", "---", "---", "---",
         "---", "---");

  {
    double zr = 50.0, zi = 0.1;
    Cpx zpz = zeta_prime_over_zeta(zr, np);
    int N_vals[] = {10, 15, 20, 25, 30, 40, 50};
    int nN = 7;
    for (int ni = 0; ni < nN; ni++) {
      int Ni = N_vals[ni];
      double l[MAXN], m[MAXN];
      for (int k = 0; k < Ni; k++)
        l[k] = ZZ[k];
      for (int k = 0; k < Ni - 1; k++)
        m[k] = 0.5 * (ZZ[k] + ZZ[k + 1]);
      double ai[MAXN], bi[MAXN];
      if (deboor(l, m, Ni, ai, bi) != 0)
        continue;
      Cpx mz = m_function(ai, bi, Ni, zr, zi);
      double delta = cpx_abs(cpx_sub(mz, zpz));
      printf("  %6d  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f\n", Ni, mz.re,
             mz.im, cpx_abs(mz), cpx_abs(zpz), delta);
    }
  }
}

static void part2(void) {
  printf("\n");
  printf("================================================================\n");
  printf("  PART 2: Krein SSF pointwise test\n");
  printf(
      "================================================================\n\n");

  int N = 30;
  double lam[MAXN], mu_arr[MAXN];
  for (int k = 0; k < N; k++)
    lam[k] = ZZ[k];
  for (int k = 0; k < N - 1; k++)
    mu_arr[k] = 0.5 * (ZZ[k] + ZZ[k + 1]);

  double a_arr[MAXN], b_arr[MAXN];
  if (deboor(lam, mu_arr, N, a_arr, b_arr) != 0) {
    printf("  dBG failed for N=%d\n", N);
    return;
  }

  int np = NPRIMES200;
  int nE = 100;
  double E_min = 10.0;
  double E_max = 140.0;
  double dE = (E_max - E_min) / (double)(nE - 1);

  printf("  xi(E) = N_zeta(E) - N_Weyl(E)\n");
  printf("  S(E) = (1/pi) * arg zeta(1/2+iE)\n");
  printf("  N_zeta = Sturm count of Jacobi eigenvalues <= E\n");
  printf("  N_Weyl = theta(E)/pi\n\n");

  printf("  %8s  %10s  %10s  %10s  %10s  %12s\n", "E", "N_zeta", "N_Weyl",
         "xi(E)", "S(E)", "|Delta|");
  printf("  %8s  %10s  %10s  %10s  %10s  %12s\n", "---", "---", "---", "---",
         "---", "---");

  for (int i = 0; i < nE; i++) {
    double E = E_min + (double)i * dE;
    int n_zeta = sturm_count(a_arr, b_arr, N, E);
    double n_weyl = theta_s(E) / M_PI;
    double xi = (double)n_zeta - n_weyl;

    double arg_z = arg_zeta_half(E, np);
    double S = arg_z / M_PI;

    double delta = fabs(xi - S);
    printf("  %8.3f  %10d  %10.4f  %10.4f  %10.4f  %12.6f\n", E, n_zeta, n_weyl,
           xi, S, delta);
  }

  printf("\n  2b: Resolvent-based xi_resolvent(E) = (1/pi) Im "
         "Tr[(J-E-ieps)^{-1}]\n");
  printf("  Using eps = 0.1, N = %d\n\n", N);

  double eps_res = 0.1;
  printf("  %8s  %10s  %10s  %12s\n", "E", "xi(E)", "xi_res(E)", "|Delta|");
  printf("  %8s  %10s  %10s  %12s\n", "---", "---", "---", "---");

  for (int i = 0; i < nE; i++) {
    double E = E_min + (double)i * dE;
    int n_zeta = sturm_count(a_arr, b_arr, N, E);
    double n_weyl = theta_s(E) / M_PI;
    double xi = (double)n_zeta - n_weyl;

    double tr_im = resolvent_trace_imag(a_arr, b_arr, N, E, eps_res);
    double xi_res = -tr_im / M_PI;

    double delta = fabs(xi - xi_res);
    printf("  %8.3f  %10.4f  %10.4f  %12.6f\n", E, xi, xi_res, delta);
  }
}

static void part3(void) {
  printf("\n");
  printf("================================================================\n");
  printf("  PART 3: Entry analytic formula — prime coeff extraction via DFT\n");
  printf(
      "================================================================\n\n");

  int N_vals[] = {30, 40, 50};
  int nN = 3;
  for (int ni = 0; ni < nN; ni++) {
    int N = N_vals[ni];
    double lam[MAXN], mu_arr[MAXN];
    for (int k = 0; k < N; k++)
      lam[k] = ZZ[k];
    for (int k = 0; k < N - 1; k++)
      mu_arr[k] = 0.5 * (ZZ[k] + ZZ[k + 1]);

    double a_arr[MAXN], b_arr[MAXN];
    if (deboor(lam, mu_arr, N, a_arr, b_arr) != 0) {
      printf("  N=%d: dBG failed\n\n", N);
      continue;
    }

    printf("  N = %d:\n", N);
    printf("  DFT of {a_k} at prime frequencies omega = log(p)\n\n");

    printf("  %6s  %10s  %12s  %12s  %12s  %12s  %10s\n", "p", "log(p)",
           "Re(c_p)", "Im(c_p)", "|c_p|", "alpha_p", "ratio");
    printf("  %6s  %10s  %12s  %12s  %12s  %12s  %10s\n", "---", "---", "---",
           "---", "---", "---", "---");

    int np_show = 30;
    if (np_show > N / 2)
      np_show = N / 2;

    double cp_list[50];
    double alpha_list[50];
    int np_used = 0;

    for (int pi = 0; pi < np_show; pi++) {
      double p = (double)PRIMES200[pi];
      double omega = log(p);
      double alpha_p = -log(p) / (2.0 * M_PI * sqrt(p));

      Cpx cp = dft_at_freq(a_arr, N, omega);
      double ratio =
          (fabs(alpha_p) > 1e-18) ? cpx_abs(cp) / fabs(alpha_p) : 0.0;

      printf("  %6.0f  %10.6f  %12.8f  %12.8f  %12.8f  %12.8f  %10.4f\n", p,
             omega, cp.re, cp.im, cpx_abs(cp), alpha_p, ratio);

      cp_list[np_used] = cpx_abs(cp);
      alpha_list[np_used] = fabs(alpha_p);
      np_used++;
    }

    double sum_xy = 0.0, sum_xx = 0.0, sum_yy = 0.0;
    double mean_x = 0.0, mean_y = 0.0;
    for (int i = 0; i < np_used; i++) {
      mean_x += cp_list[i];
      mean_y += alpha_list[i];
    }
    mean_x /= (double)np_used;
    mean_y /= (double)np_used;
    for (int i = 0; i < np_used; i++) {
      double dx = cp_list[i] - mean_x;
      double dy = alpha_list[i] - mean_y;
      sum_xy += dx * dy;
      sum_xx += dx * dx;
      sum_yy += dy * dy;
    }
    double corr =
        (sum_xx * sum_yy > 1e-30) ? sum_xy / sqrt(sum_xx * sum_yy) : 0.0;

    double lambda_num = 0.0, lambda_den = 0.0;
    for (int i = 0; i < np_used; i++) {
      lambda_num += cp_list[i] * alpha_list[i];
      lambda_den += alpha_list[i] * alpha_list[i];
    }
    double lambda_opt = (lambda_den > 1e-30) ? lambda_num / lambda_den : 0.0;

    double ss_res = 0.0, ss_tot = 0.0;
    for (int i = 0; i < np_used; i++) {
      double pred = lambda_opt * alpha_list[i];
      double err = cp_list[i] - pred;
      ss_res += err * err;
      ss_tot += (cp_list[i] - mean_x) * (cp_list[i] - mean_x);
    }
    double r2 = (ss_tot > 1e-30) ? 1.0 - ss_res / ss_tot : 0.0;

    printf("\n  LSQ fit: c_p ~ lambda * alpha_p\n");
    printf("  lambda_opt = %.6f\n", lambda_opt);
    printf("  Pearson r = %.6f\n", corr);
    printf("  R^2 = %.6f\n\n", r2);

    printf("  Reconstruction R^2: a_k = Sum_p [Re(c_p)*cos(k ln p) - "
           "Im(c_p)*sin(k ln p)]\n");
    double ss_res2 = 0.0, ss_tot2 = 0.0;
    double mean_a = 0.0;
    for (int k = 0; k < N; k++)
      mean_a += a_arr[k];
    mean_a /= (double)N;
    for (int k = 0; k < N; k++) {
      double pred = 0.0;
      for (int pi = 0; pi < np_show; pi++) {
        double p = (double)PRIMES200[pi];
        double omega = log(p);
        Cpx cp = dft_at_freq(a_arr, N, omega);
        pred += cp.re * cos(omega * (double)k) - cp.im * sin(omega * (double)k);
      }
      pred *= (double)N;
      double err = a_arr[k] - pred;
      ss_res2 += err * err;
      ss_tot2 += (a_arr[k] - mean_a) * (a_arr[k] - mean_a);
    }
    double r2_recon = (ss_tot2 > 1e-30) ? 1.0 - ss_res2 / ss_tot2 : 0.0;
    printf("  R^2 (DFT reconstruction of a_k) = %.6f\n", r2_recon);

    printf("\n  DFT of {b_k} at prime frequencies:\n\n");
    printf("  %6s  %10s  %12s  %12s  %12s\n", "p", "log(p)", "Re(c_p)",
           "Im(c_p)", "|c_p|");
    printf("  %6s  %10s  %12s  %12s  %12s\n", "---", "---", "---", "---",
           "---");

    for (int pi = 0; pi < np_show; pi++) {
      double p = (double)PRIMES200[pi];
      double omega = log(p);
      Cpx cp = dft_at_freq(b_arr, N - 1, omega);
      printf("  %6.0f  %10.6f  %12.8f  %12.8f  %12.8f\n", p, omega, cp.re,
             cp.im, cpx_abs(cp));
    }

    double ss_res_b = 0.0, ss_tot_b = 0.0;
    double mean_b = 0.0;
    for (int k = 0; k < N - 1; k++)
      mean_b += b_arr[k];
    mean_b /= (double)(N - 1);
    for (int k = 0; k < N - 1; k++) {
      double pred = 0.0;
      for (int pi = 0; pi < np_show; pi++) {
        double p = (double)PRIMES200[pi];
        double omega = log(p);
        Cpx cp = dft_at_freq(b_arr, N - 1, omega);
        pred += cp.re * cos(omega * (double)k) - cp.im * sin(omega * (double)k);
      }
      pred *= (double)(N - 1);
      double err = b_arr[k] - pred;
      ss_res_b += err * err;
      ss_tot_b += (b_arr[k] - mean_b) * (b_arr[k] - mean_b);
    }
    double r2_b = (ss_tot_b > 1e-30) ? 1.0 - ss_res_b / ss_tot_b : 0.0;
    printf("\n  R^2 (DFT reconstruction of b_k) = %.6f\n\n", r2_b);

    printf("  Ratio c_p/alpha_p constancy check:\n");
    printf("  %6s  %12s  %12s  %12s\n", "p", "|c_p|", "alpha_p", "ratio");
    printf("  %6s  %12s  %12s  %12s\n", "---", "---", "---", "---");
    for (int pi = 0; pi < np_show; pi++) {
      double p = (double)PRIMES200[pi];
      double omega = log(p);
      double alpha_p = -log(p) / (2.0 * M_PI * sqrt(p));
      Cpx cp = dft_at_freq(a_arr, N, omega);
      double ratio =
          (fabs(alpha_p) > 1e-18) ? cpx_abs(cp) / fabs(alpha_p) : 0.0;
      printf("  %6.0f  %12.8f  %12.8f  %12.4f\n", p, cpx_abs(cp), fabs(alpha_p),
             ratio);
    }
    printf("\n");
  }
}

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

  printf("================================================================\n");
  printf("  m-function Bridge: Three Critical Computations for RH Proof\n");
  printf("  Part 1: m(z) vs zeta'/zeta direct comparison\n");
  printf("  Part 2: Krein SSF pointwise test\n");
  printf("  Part 3: Entry analytic formula — prime coeff extraction\n");
  printf("  Using 200 primes, Euler exponential form\n");
  printf("================================================================\n");

  part1();
  part2();
  part3();

  printf("\nDone.\n");
  return 0;
}
