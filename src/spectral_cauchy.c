/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Spectral measure Cauchy transform
 * @paper   yamaguchi-rh-2026.tex, §10.6
 * @theorem Theorem III
 * @proof   G(z) = sum w_k/(gamma_k-z) vs zeta'/zeta
 * @step    4
 *
 * spectral_cauchy.c — Spectral measure Cauchy transform G(z) vs zeta'/zeta
 *
 * KEY INSIGHT: m(z) = <e1|(J-z)^{-1}|e1> is just the (1,1) resolvent entry,
 * which converges to 0 as N->inf. The CORRECT object is the full spectral
 * measure Cauchy transform G(z) = sum_k w_k/(gamma_k - z) where w_k are
 * the dBG spectral weights.
 *
 * Part 1: G(z) at zeta zeros — residue = w_k
 * Part 2: G(z) vs zeta'/zeta at regular points — convergence with N
 * Part 3: Weighted vs Uniform vs m(z) — which matches zeta'/zeta?
 * Part 4: Residue check — w_k vs 1/|zeta'(rho_k)|
 *
 * Compile: gcc -Wall -Wextra -Wconversion -Wshadow -Werror -O3
 *          -fno-strict-aliasing -fno-peel-loops -fno-unswitch-loops
 *          -Isrc -o spectral_cauchy.exe src/spectral_cauchy.c -lm
 */

#include "random_matrix_utils.h"
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
#define MAXN 2100
#define NPRIMES 200

static const int PRIMES[NPRIMES] = {
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
                  double *b, double *w_out) {
  for (int k = 0; k < N - 1; k++)
    if (!(lam[k] < mu[k] && mu[k] < lam[k + 1]))
      return -1;
  double w[MAXN];
  double ws = 0.0;
  for (int k = 0; k < N; k++) {
    double log_nm = 0.0;
    double log_dn = 0.0;
    int sign_nm = 1;
    int sign_dn = 1;
    for (int j = 0; j < N - 1; j++) {
      double d = lam[k] - mu[j];
      if (d == 0.0)
        return -1;
      if (d < 0.0) {
        sign_nm = -sign_nm;
        d = -d;
      }
      log_nm += log(d);
    }
    for (int j = 0; j < N; j++) {
      if (j == k)
        continue;
      double d = lam[k] - lam[j];
      if (d == 0.0)
        return -1;
      if (d < 0.0) {
        sign_dn = -sign_dn;
        d = -d;
      }
      log_dn += log(d);
    }
    double log_w = log_nm - log_dn;
    if (log_w > 700.0)
      return -2;
    w[k] = (double)(sign_nm * sign_dn) * exp(log_w);
    if (w[k] < 0.0)
      return -2;
    ws += w[k];
  }
  for (int k = 0; k < N; k++)
    w[k] /= ws;
  if (w_out) {
    for (int k = 0; k < N; k++)
      w_out[k] = w[k];
  }
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

static Cpx zeta_prime_over_zeta(double E, int np) {
  Cpx sum = cpx_make(0.0, 0.0);
  for (int pi = 0; pi < np; pi++) {
    double p = (double)PRIMES[pi];
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

static Cpx zeta_approx(double E, int np) {
  Cpx prod = cpx_make(1.0, 0.0);
  for (int pi = 0; pi < np; pi++) {
    double p = (double)PRIMES[pi];
    Cpx phase = cpx_p_ik(p, E);
    Cpx ps = cpx_scale(sqrt(p), phase);
    Cpx factor = cpx_sub(cpx_make(1.0, 0.0), cpx_inv(ps));
    prod = cpx_mul(prod, factor);
  }
  return prod;
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

static Cpx spectral_cauchy(const double *lam, const double *w, int N, double zr,
                           double zi) {
  Cpx z = cpx_make(zr, zi);
  Cpx sum = cpx_make(0.0, 0.0);
  for (int k = 0; k < N; k++) {
    Cpx denom = cpx_sub(cpx_make(lam[k], 0.0), z);
    if (cpx_abs(denom) < 1e-14) {
      sum = cpx_add(sum, cpx_make(w[k] * 1e14, 0.0));
    } else {
      sum = cpx_add(sum, cpx_scale(w[k], cpx_inv(denom)));
    }
  }
  return sum;
}

static Cpx stieltjes_uniform(const double *lam, int N, double zr, double zi) {
  Cpx z = cpx_make(zr, zi);
  Cpx sum = cpx_make(0.0, 0.0);
  for (int k = 0; k < N; k++) {
    Cpx denom = cpx_sub(cpx_make(lam[k], 0.0), z);
    if (cpx_abs(denom) < 1e-14)
      continue;
    sum = cpx_add(sum, cpx_inv(denom));
  }
  return cpx_scale(1.0 / (double)N, sum);
}

__attribute__((unused)) static double theta_s(double t) {
  if (t <= 1.0)
    return -M_PI / 8.0;
  double x = t / (2.0 * M_PI);
  double u = 1.0 / t;
  return 0.5 * t * log(x) - 0.5 * t - M_PI / 8.0 + u / 48.0 +
         7.0 * u * u * u / 5760.0;
}

static void part1(void) {
  printf("\n");
  printf("================================================================\n");
  // G(z) -> 0, not zeta'/z. Paper Section 8.5, Theorem III.
  printf("  PART 1: Spectral Cauchy G(z) at zeta zeros — residue check\n");
  printf(
      "================================================================\n\n");

  printf("  G(z) = sum_k w_k/(gamma_k - z),  z = gamma_j + i*eps\n");
  printf("  At a zero, residue(G, gamma_j) = w_j\n\n");

  int N_vals[] = {10, 20, 30, 50, 100, 200};
  int nN = 6;
  double eps_vals[] = {0.1, 0.01, 0.001};
  int n_eps = 3;

  for (int ni = 0; ni < nN; ni++) {
    int N = N_vals[ni];
    double lam[MAXN], mu_arr[MAXN];
    for (int k = 0; k < N; k++)
      lam[k] = ZETA_ZEROS_2000[k];
    for (int k = 0; k < N - 1; k++)
      mu_arr[k] = 0.5 * (ZETA_ZEROS_2000[k] + ZETA_ZEROS_2000[k + 1]);

    double a_arr[MAXN], b_arr[MAXN], w[MAXN];
    if (deboor(lam, mu_arr, N, a_arr, b_arr, w) != 0) {
      printf("  N=%d: dBG failed\n\n", N);
      continue;
    }

    printf("  N = %d:\n", N);
    printf("  %4s  %10s  %5s  %12s  %12s  %12s  %12s  %12s\n", "j", "E", "eps",
           "|G(z)|", "w_j", "|z'/z|", "Res_check", "ratio");
    printf("  %4s  %10s  %5s  %12s  %12s  %12s  %12s  %12s\n", "---", "---",
           "---", "---", "---", "---", "---", "---");

    int j_vals[] = {0, 1, 2, 4, 9};
    int nj = 5;
    if (N <= 10) {
      j_vals[3] = 3;
      j_vals[4] = 4;
      nj = 5;
    }

    for (int ji = 0; ji < nj; ji++) {
      int j = j_vals[ji];
      if (j >= N)
        continue;
      double E = ZETA_ZEROS_2000[j];
      for (int ei = 0; ei < n_eps; ei++) {
        double eps = eps_vals[ei];
        Cpx gz = spectral_cauchy(lam, w, N, E, eps);
        Cpx zpz = zeta_prime_over_zeta(E, NPRIMES);
        double res_check = cpx_abs(gz) * eps;
        double ratio =
            (cpx_abs(zpz) > 1e-10) ? cpx_abs(gz) / cpx_abs(zpz) : 0.0;
        printf("  %4d  %10.4f  %5.3f  %12.4f  %12.6f  %12.4f  %12.6f  %12.4f\n",
               j, E, eps, cpx_abs(gz), w[j], cpx_abs(zpz), res_check, ratio);
      }
    }
    printf("\n");
  }
}

static void part2(void) {
  printf("\n");
  printf("================================================================\n");
  printf("  PART 2: G(z) vs zeta'/zeta at regular points — convergence\n");
  printf(
      "================================================================\n\n");

  printf("  G(z) = sum_k w_k/(gamma_k - z)  vs  zeta'/zeta(1/2+iz)\n");
  printf("  At regular points (midpoints between zeros)\n\n");

  int np_vals[] = {50, 100, 200};
  int n_np = 3;

  for (int npi = 0; npi < n_np; npi++) {
    int np_use = np_vals[npi];
    printf("  Using %d primes for Euler product:\n", np_use);
    printf("  %6s  %12s  %12s  %12s  %12s  %12s\n", "N", "|G(z)|", "|z'/z|",
           "ratio", "Re(G)", "Im(G)");
    printf("  %6s  %12s  %12s  %12s  %12s  %12s\n", "---", "---", "---", "---",
           "---", "---");

    int N_vals[] = {10, 20, 30, 50, 100, 200, 500, 1000};
    int nN = 8;

    for (int ni = 0; ni < nN; ni++) {
      int N = N_vals[ni];
      double lam[MAXN], mu_arr[MAXN];
      for (int k = 0; k < N; k++)
        lam[k] = ZETA_ZEROS_2000[k];
      for (int k = 0; k < N - 1; k++)
        mu_arr[k] = 0.5 * (ZETA_ZEROS_2000[k] + ZETA_ZEROS_2000[k + 1]);

      double a_arr[MAXN], b_arr[MAXN], w[MAXN];
      if (deboor(lam, mu_arr, N, a_arr, b_arr, w) != 0)
        continue;

      double zr = 0.5 * (ZETA_ZEROS_2000[9] + ZETA_ZEROS_2000[10]);
      double zi = 0.1;

      Cpx gz = spectral_cauchy(lam, w, N, zr, zi);
      Cpx zpz = zeta_prime_over_zeta(zr, np_use);
      double ratio = (cpx_abs(zpz) > 1e-10) ? cpx_abs(gz) / cpx_abs(zpz) : 0.0;

      printf("  %6d  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f\n", N, cpx_abs(gz),
             cpx_abs(zpz), ratio, gz.re, gz.im);
    }
    printf("\n");
  }

  printf(
      "  Convergence of ratio G(z)/zeta'/zeta at multiple regular points:\n");
  printf("  %6s  %10s  %12s  %12s  %12s  %12s\n", "N", "E_mid", "|G(z)|",
         "|z'/z|", "ratio", "|G-S_unif|/|G|");
  printf("  %6s  %10s  %12s  %12s  %12s  %12s\n", "---", "---", "---", "---",
         "---", "---");

  int mid_indices[] = {4, 9, 14, 19, 29, 49};
  int n_mid = 6;
  int N_test = 200;

  double lam[MAXN], mu_arr[MAXN];
  for (int k = 0; k < N_test; k++)
    lam[k] = ZETA_ZEROS_2000[k];
  for (int k = 0; k < N_test - 1; k++)
    mu_arr[k] = 0.5 * (ZETA_ZEROS_2000[k] + ZETA_ZEROS_2000[k + 1]);

  double a_arr[MAXN], b_arr[MAXN], w[MAXN];
  if (deboor(lam, mu_arr, N_test, a_arr, b_arr, w) == 0) {
    for (int mi = 0; mi < n_mid; mi++) {
      int idx = mid_indices[mi];
      if (idx >= N_test - 1)
        continue;
      double zr = 0.5 * (ZETA_ZEROS_2000[idx] + ZETA_ZEROS_2000[idx + 1]);
      double zi = 0.1;

      Cpx gz = spectral_cauchy(lam, w, N_test, zr, zi);
      Cpx su = stieltjes_uniform(lam, N_test, zr, zi);
      Cpx zpz = zeta_prime_over_zeta(zr, NPRIMES);
      double ratio = (cpx_abs(zpz) > 1e-10) ? cpx_abs(gz) / cpx_abs(zpz) : 0.0;
      double rel_diff =
          (cpx_abs(gz) > 1e-10) ? cpx_abs(cpx_sub(gz, su)) / cpx_abs(gz) : 0.0;

      printf("  %6d  %10.4f  %12.6f  %12.6f  %12.6f  %12.6f\n", N_test, zr,
             cpx_abs(gz), cpx_abs(zpz), ratio, rel_diff);
    }
  }
  printf("\n");
}

static void part3(void) {
  printf("\n");
  printf("================================================================\n");
  printf("  PART 3: Weighted vs Uniform vs m(z) — which matches zeta'/zeta?\n");
  printf(
      "================================================================\n\n");

  double zr_vals[] = {50.0, 100.0, 150.0};
  double zi = 0.1;
  int nz = 3;

  printf("  z = E + 0.1i\n\n");

  for (int zi_idx = 0; zi_idx < nz; zi_idx++) {
    double zr = zr_vals[zi_idx];
    printf("  z = %.1f + 0.1i:\n", zr);
    printf("  %6s  %12s  %12s  %12s  %12s  %12s  %12s\n", "N", "|G_weight|",
           "|S_uniform|", "|m(z)|", "|z'/z|", "rG", "rS");
    printf("  %6s  %12s  %12s  %12s  %12s  %12s  %12s\n", "---", "---", "---",
           "---", "---", "---", "---");

    Cpx zpz = zeta_prime_over_zeta(zr, NPRIMES);

    int N_vals[] = {10, 20, 30, 50, 100, 200};
    int nN = 6;

    for (int ni = 0; ni < nN; ni++) {
      int N = N_vals[ni];
      double lam[MAXN], mu_arr[MAXN];
      for (int k = 0; k < N; k++)
        lam[k] = ZETA_ZEROS_2000[k];
      for (int k = 0; k < N - 1; k++)
        mu_arr[k] = 0.5 * (ZETA_ZEROS_2000[k] + ZETA_ZEROS_2000[k + 1]);

      double a_arr[MAXN], b_arr[MAXN], w[MAXN];
      if (deboor(lam, mu_arr, N, a_arr, b_arr, w) != 0)
        continue;

      Cpx gz = spectral_cauchy(lam, w, N, zr, zi);
      Cpx su = stieltjes_uniform(lam, N, zr, zi);
      Cpx mz = m_function(a_arr, b_arr, N, zr, zi);

      double rG = (cpx_abs(zpz) > 1e-10) ? cpx_abs(gz) / cpx_abs(zpz) : 0.0;
      double rS = (cpx_abs(zpz) > 1e-10) ? cpx_abs(su) / cpx_abs(zpz) : 0.0;

      printf("  %6d  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f\n", N,
             cpx_abs(gz), cpx_abs(su), cpx_abs(mz), cpx_abs(zpz), rG, rS);
    }
    printf("\n");
  }

  printf("  Full complex comparison at z = 50.0 + 0.1i, N=50:\n");
  printf("  %12s  %12s  %12s  %12s  %12s  %12s\n", "", "Re", "Im", "|.|",
         "arg(deg)", "");
  printf("  %12s  %12s  %12s  %12s  %12s  %12s\n", "---", "---", "---", "---",
         "---", "---");
  {
    int N = 50;
    double lam[MAXN], mu_arr[MAXN];
    for (int k = 0; k < N; k++)
      lam[k] = ZETA_ZEROS_2000[k];
    for (int k = 0; k < N - 1; k++)
      mu_arr[k] = 0.5 * (ZETA_ZEROS_2000[k] + ZETA_ZEROS_2000[k + 1]);

    double a_arr[MAXN], b_arr[MAXN], w[MAXN];
    if (deboor(lam, mu_arr, N, a_arr, b_arr, w) == 0) {
      double zr = 50.0;
      Cpx gz = spectral_cauchy(lam, w, N, zr, zi);
      Cpx su = stieltjes_uniform(lam, N, zr, zi);
      Cpx mz = m_function(a_arr, b_arr, N, zr, zi);
      Cpx zpz = zeta_prime_over_zeta(zr, NPRIMES);

      printf("  %12s  %12.6f  %12.6f  %12.6f  %12.4f\n", "G_weighted", gz.re,
             gz.im, cpx_abs(gz), cpx_arg(gz) * 180.0 / M_PI);
      printf("  %12s  %12.6f  %12.6f  %12.6f  %12.4f\n", "S_uniform", su.re,
             su.im, cpx_abs(su), cpx_arg(su) * 180.0 / M_PI);
      printf("  %12s  %12.6f  %12.6f  %12.6f  %12.4f\n", "m(z)", mz.re, mz.im,
             cpx_abs(mz), cpx_arg(mz) * 180.0 / M_PI);
      printf("  %12s  %12.6f  %12.6f  %12.6f  %12.4f\n", "zeta'/zeta", zpz.re,
             zpz.im, cpx_abs(zpz), cpx_arg(zpz) * 180.0 / M_PI);
    }
  }
  printf("\n");
}

static void part4(void) {
  printf("\n");
  printf("================================================================\n");
  printf("  PART 4: Residue check — w_k vs 1/|zeta'(rho_k)|\n");
  printf(
      "================================================================\n\n");

  printf("  If the HP operator is correct, the spectral measure dmu = sum w_k "
         "delta(gamma_k)\n");
  printf("  should relate to the zeros of zeta via w_k = C / |zeta'(rho_k)| "
         "for some C.\n\n");

  int N_vals[] = {10, 20, 30, 50, 100};
  int nN = 5;

  for (int ni = 0; ni < nN; ni++) {
    int N = N_vals[ni];
    double lam[MAXN], mu_arr[MAXN];
    for (int k = 0; k < N; k++)
      lam[k] = ZETA_ZEROS_2000[k];
    for (int k = 0; k < N - 1; k++)
      mu_arr[k] = 0.5 * (ZETA_ZEROS_2000[k] + ZETA_ZEROS_2000[k + 1]);

    double a_arr[MAXN], b_arr[MAXN], w[MAXN];
    if (deboor(lam, mu_arr, N, a_arr, b_arr, w) != 0) {
      printf("  N=%d: dBG failed\n\n", N);
      continue;
    }

    printf("  N = %d:\n", N);
    printf("  %4s  %12s  %12s  %12s  %12s  %12s\n", "k", "gamma_k", "w_k",
           "|zeta'(eps)|", "1/|z'|*C", "ratio w/(1/|z'|)");
    printf("  %4s  %12s  %12s  %12s  %12s  %12s\n", "---", "---", "---", "---",
           "---", "---");

    double eps_fd = 0.01;
    double inv_zprime[MAXN];
    double sum_inv = 0.0;
    double sum_w = 0.0;

    for (int k = 0; k < N; k++) {
      Cpx zp = zeta_approx(ZETA_ZEROS_2000[k] + eps_fd, NPRIMES);
      Cpx zm = zeta_approx(ZETA_ZEROS_2000[k] - eps_fd, NPRIMES);
      Cpx deriv = cpx_scale(1.0 / (2.0 * eps_fd), cpx_sub(zp, zm));
      inv_zprime[k] = 1.0 / cpx_abs(deriv);
      sum_inv += inv_zprime[k];
      sum_w += w[k];
    }

    double C_norm = (sum_inv > 1e-30) ? sum_w / sum_inv : 1.0;

    int show_n = (N <= 20) ? N : 15;
    for (int k = 0; k < show_n; k++) {
      double scaled = inv_zprime[k] * C_norm;
      double ratio = (scaled > 1e-18) ? w[k] / scaled : 0.0;
      printf("  %4d  %12.4f  %12.6f  %12.6f  %12.6f  %12.6f\n", k,
             ZETA_ZEROS_2000[k], w[k], 1.0 / inv_zprime[k], scaled, ratio);
    }

    double ss_res = 0.0, ss_tot = 0.0;
    for (int k = 0; k < N; k++) {
      double pred = inv_zprime[k] * C_norm;
      double err = w[k] - pred;
      ss_res += err * err;
      ss_tot += (w[k] - sum_w / (double)N) * (w[k] - sum_w / (double)N);
    }
    double r2 = (ss_tot > 1e-30) ? 1.0 - ss_res / ss_tot : 0.0;

    printf("  Normalization C = %.6f (so that sum(w_k) = C * sum(1/|z'|))\n",
           C_norm);
    printf("  R^2 (w_k vs C/|zeta'(rho_k)|) = %.6f\n", r2);
    printf("  sum(w_k) = %.6f,  C*sum(1/|z'|) = %.6f\n\n", sum_w,
           C_norm * sum_inv);
  }

  printf("  Convergence of R^2 with N:\n");
  printf("  %6s  %12s  %12s  %12s\n", "N", "R^2", "C", "sum(w)/N");
  printf("  %6s  %12s  %12s  %12s\n", "---", "---", "---", "---");

  int N_conv[] = {5, 8, 10, 15, 20, 30, 50, 100, 200, 500, 1000};
  int nNc = 11;

  for (int ni = 0; ni < nNc; ni++) {
    int N = N_conv[ni];
    double lam[MAXN], mu_arr[MAXN];
    for (int k = 0; k < N; k++)
      lam[k] = ZETA_ZEROS_2000[k];
    for (int k = 0; k < N - 1; k++)
      mu_arr[k] = 0.5 * (ZETA_ZEROS_2000[k] + ZETA_ZEROS_2000[k + 1]);

    double a_arr[MAXN], b_arr[MAXN], w[MAXN];
    if (deboor(lam, mu_arr, N, a_arr, b_arr, w) != 0)
      continue;

    double eps_fd2 = 0.01;
    double inv_zp[MAXN];
    double si = 0.0, sw = 0.0;
    for (int k = 0; k < N; k++) {
      Cpx zp = zeta_approx(ZETA_ZEROS_2000[k] + eps_fd2, NPRIMES);
      Cpx zm = zeta_approx(ZETA_ZEROS_2000[k] - eps_fd2, NPRIMES);
      Cpx deriv = cpx_scale(1.0 / (2.0 * eps_fd2), cpx_sub(zp, zm));
      inv_zp[k] = 1.0 / cpx_abs(deriv);
      si += inv_zp[k];
      sw += w[k];
    }
    double C = (si > 1e-30) ? sw / si : 1.0;

    double sr = 0.0, st = 0.0;
    for (int k = 0; k < N; k++) {
      double pred = inv_zp[k] * C;
      double err = w[k] - pred;
      sr += err * err;
      st += (w[k] - sw / (double)N) * (w[k] - sw / (double)N);
    }
    double r2v = (st > 1e-30) ? 1.0 - sr / st : 0.0;

    printf("  %6d  %12.6f  %12.6f  %12.6f\n", N, r2v, C, sw / (double)N);
  }
}

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

  printf("================================================================\n");
  // G(z) -> 0, not zeta'/z. Paper Section 8.5, Theorem III.
  printf("  Spectral Cauchy Transform: G(z) = sum_k w_k/(gamma_k - z)\n");
  printf("  vs zeta'/zeta(1/2+iz)\n");
  printf("  Using %d zeta zeros, %d primes, Euler exponential form\n",
         N_REF_2000, NPRIMES);
  printf("  Part 1: Residue check at zeta zeros\n");
  printf("  Part 2: G(z) vs zeta'/zeta at regular points\n");
  printf("  Part 3: Weighted vs Uniform vs m(z)\n");
  printf("  Part 4: w_k vs 1/|zeta'(rho_k)|\n");
  printf("================================================================\n");

  part1();
  part2();
  part3();
  part4();

  printf("\nDone.\n");
  return 0;
}
