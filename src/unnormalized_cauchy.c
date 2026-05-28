/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Unnormalized Cauchy transform
 * @paper   yamaguchi-rh-2026.tex, §10.6
 * @theorem Theorem III
 * @proof   Obstruction analysis
 * @step    4
 *
 * unnormalized_cauchy.c â€” Raw (unnormalized) Lagrange weight Cauchy transform
 *
 * Tests whether G_raw(z) = sum_k w_k^raw/(gamma_k - z) converges to
 * zeta'/zeta(z) as N -> inf, where w_k^raw are the unnormalized dBG spectral
 * weights.
 *
 * Computations:
 *   1. Raw weight magnitude analysis
 *   2. G_raw vs zeta'/zeta at regular points
 *   3. Normalized vs Raw near zeta zeros
 *   4. Convergence of ratio G_raw/zeta'/zeta across multiple points
 *
 * Compile: gcc -Wall -Wextra -Wconversion -Wshadow -Werror -O3
 *          -fno-strict-aliasing -fno-peel-loops -fno-unswitch-loops
 *          -Isrc -o unnormalized_cauchy.exe src/unnormalized_cauchy.c -lm
 */

#include "random_matrix_utils.h"
#include "refdata_2000.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

static double logsumexp(const double *vals, int n) {
  if (n <= 0)
    return -1e300;
  double mx = vals[0];
  for (int i = 1; i < n; i++) {
    if (vals[i] > mx)
      mx = vals[i];
  }
  double s = 0.0;
  for (int i = 0; i < n; i++) {
    s += exp(vals[i] - mx);
  }
  return mx + log(s);
}

static void compute_raw_weights(const double *lam, const double *mu, int N,
                                double *log_w, int *sign_w) {
  for (int k = 0; k < N; k++) {
    double log_num = 0.0;
    double log_den = 0.0;
    int s_num = 1;
    int s_den = 1;
    for (int j = 0; j < N - 1; j++) {
      double d = lam[k] - mu[j];
      if (d < 0.0) {
        s_num = -s_num;
        d = -d;
      }
      log_num += log(d);
    }
    for (int j = 0; j < N; j++) {
      if (j == k)
        continue;
      double d = lam[k] - lam[j];
      if (d < 0.0) {
        s_den = -s_den;
        d = -d;
      }
      log_den += log(d);
    }
    log_w[k] = log_num - log_den;
    sign_w[k] = s_num * s_den;
  }
}

static Cpx raw_cauchy_log(const double *lam, const double *log_w,
                          const int *sign_w, int N, double zr, double zi) {
  double lpos_re[MAXN], lneg_re[MAXN];
  double lpos_im[MAXN], lneg_im[MAXN];
  int npr = 0, nnr = 0, npi = 0, nni = 0;
  double log_abs_zi = log(fabs(zi));
  int sgn_zi = (zi >= 0.0) ? 1 : -1;

  for (int k = 0; k < N; k++) {
    double dr = lam[k] - zr;
    double dist2 = dr * dr + zi * zi;
    double log_dist2 = log(dist2);

    if (fabs(dr) > 1e-30) {
      double lm = log_w[k] + log(fabs(dr)) - log_dist2;
      int sr = sign_w[k] * ((dr > 0.0) ? 1 : -1);
      if (sr > 0)
        lpos_re[npr++] = lm;
      else
        lneg_re[nnr++] = lm;
    }

    {
      double lm = log_w[k] + log_abs_zi - log_dist2;
      int si = sign_w[k] * sgn_zi;
      if (si > 0)
        lpos_im[npi++] = lm;
      else
        lneg_im[nni++] = lm;
    }
  }

  double pr_v = logsumexp(lpos_re, npr);
  double nr_v = logsumexp(lneg_re, nnr);
  double pi_v = logsumexp(lpos_im, npi);
  double ni_v = logsumexp(lneg_im, nni);

  double re;
  {
    double mx = (pr_v > nr_v) ? pr_v : nr_v;
    if (pr_v < -1e290 && nr_v < -1e290) {
      re = 0.0;
    } else if (nr_v < -1e290) {
      re = exp(pr_v);
    } else if (pr_v < -1e290) {
      re = -exp(nr_v);
    } else {
      re = exp(mx) * (exp(pr_v - mx) - exp(nr_v - mx));
    }
  }
  double im;
  {
    double mx = (pi_v > ni_v) ? pi_v : ni_v;
    if (pi_v < -1e290 && ni_v < -1e290) {
      im = 0.0;
    } else if (ni_v < -1e290) {
      im = exp(pi_v);
    } else if (pi_v < -1e290) {
      im = -exp(ni_v);
    } else {
      im = exp(mx) * (exp(pi_v - mx) - exp(ni_v - mx));
    }
  }

  return cpx_make(re, im);
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

static Cpx norm_cauchy(const double *lam, const double *w, int N, double zr,
                       double zi) {
  Cpx z = cpx_make(zr, zi);
  Cpx s = cpx_make(0.0, 0.0);
  for (int k = 0; k < N; k++) {
    Cpx denom = cpx_sub(cpx_make(lam[k], 0.0), z);
    if (cpx_abs(denom) < 1e-14) {
      s = cpx_add(s, cpx_make(w[k] * 1e14, 0.0));
    } else {
      s = cpx_add(s, cpx_scale(w[k], cpx_inv(denom)));
    }
  }
  return s;
}

static void computation1(void) {
  printf("\n");
  printf("================================================================\n");
  printf("  COMPUTATION 1: Raw weights analysis\n");
  printf(
      "================================================================\n\n");

  int N_vals[] = {10, 20, 30, 50, 100, 200, 500, 1000};
  int nN = 8;

  for (int ni = 0; ni < nN; ni++) {
    int N = N_vals[ni];
    double lam[MAXN], mu_arr[MAXN];
    for (int k = 0; k < N; k++)
      lam[k] = ZETA_ZEROS_2000[k];
    for (int k = 0; k < N - 1; k++)
      mu_arr[k] = 0.5 * (ZETA_ZEROS_2000[k] + ZETA_ZEROS_2000[k + 1]);

    double log_w[MAXN];
    int sign_w[MAXN];
    compute_raw_weights(lam, mu_arr, N, log_w, sign_w);

    double min_lw = log_w[0], max_lw = log_w[0];
    for (int k = 1; k < N; k++) {
      if (log_w[k] < min_lw)
        min_lw = log_w[k];
      if (log_w[k] > max_lw)
        max_lw = log_w[k];
    }
    double log_sum = logsumexp(log_w, N);

    printf("  N = %d:\n", N);
    if (log_sum < 700.0) {
      printf("    sum(w_raw) = %.6e  (log = %.4f)\n", exp(log_sum), log_sum);
    } else {
      printf("    sum(w_raw) = OVERFLOW  (log = %.4f)\n", log_sum);
    }
    printf("    log10|w_raw| range: [%.4f, %.4f]  |max/min| = %.4e\n",
           min_lw / log(10.0), max_lw / log(10.0), exp(max_lw - min_lw));

    int ri0 = 0, ri1 = N / 4, ri2 = N / 2, ri3 = 3 * N / 4, ri4 = N - 1;
    printf(
        "    k=%d: log10|w|=%.4f   k=%d: log10|w|=%.4f   k=%d: log10|w|=%.4f\n",
        ri0, log_w[ri0] / log(10.0), ri1, log_w[ri1] / log(10.0), ri2,
        log_w[ri2] / log(10.0));
    printf("    k=%d: log10|w|=%.4f   k=%d: log10|w|=%.4f\n", ri3,
           log_w[ri3] / log(10.0), ri4, log_w[ri4] / log(10.0));
    printf("\n");
  }
}

static void computation2(void) {
  printf("\n");
  printf("================================================================\n");
  printf("  COMPUTATION 2: Raw Cauchy transform vs zeta'/zeta at regular "
         "points\n");
  printf(
      "================================================================\n\n");

  int N_vals[] = {10, 20, 30, 50, 100, 200, 500, 1000};
  int nN = 8;
  double zi = 0.1;

  double zr_a = 0.5 * (ZETA_ZEROS_2000[9] + ZETA_ZEROS_2000[10]);
  double zr_b = 100.0;
  double zr_c = 35.0;
  double zr_pts[3];
  zr_pts[0] = zr_a;
  zr_pts[1] = zr_b;
  zr_pts[2] = zr_c;

  for (int pi = 0; pi < 3; pi++) {
    double zr = zr_pts[pi];
    Cpx zpz = zeta_prime_over_zeta(zr, NPRIMES);
    printf("  z = %.4f + %.1fi   |z'/z| = %.6f\n", zr, zi, cpx_abs(zpz));
    printf("  %6s  %12s  %12s  %12s  %12s  %12s\n", "N", "|G_raw|", "|z'/z|",
           "ratio", "Re(G)", "Im(G)");
    printf("  %6s  %12s  %12s  %12s  %12s  %12s\n", "---", "---", "---", "---",
           "---", "---");

    for (int ni = 0; ni < nN; ni++) {
      int N = N_vals[ni];
      double lam[MAXN], mu_arr[MAXN];
      for (int k = 0; k < N; k++)
        lam[k] = ZETA_ZEROS_2000[k];
      for (int k = 0; k < N - 1; k++)
        mu_arr[k] = 0.5 * (ZETA_ZEROS_2000[k] + ZETA_ZEROS_2000[k + 1]);

      double log_w[MAXN];
      int sign_w[MAXN];
      compute_raw_weights(lam, mu_arr, N, log_w, sign_w);

      Cpx gr = raw_cauchy_log(lam, log_w, sign_w, N, zr, zi);
      double ratio = (cpx_abs(zpz) > 1e-15) ? cpx_abs(gr) / cpx_abs(zpz) : 0.0;
      printf("  %6d  %12.4e  %12.6f  %12.6f  %12.4e  %12.4e\n", N, cpx_abs(gr),
             cpx_abs(zpz), ratio, gr.re, gr.im);
    }
    printf("\n");
  }
}

static void computation3(void) {
  printf("\n");
  printf("================================================================\n");
  printf("  COMPUTATION 3: Normalized vs Raw at zeta zeros\n");
  printf(
      "================================================================\n\n");

  int j_vals[] = {0, 4, 9};
  int nj = 3;
  double eps_vals[] = {0.1, 0.01};
  int neps = 2;

  int N_vals[] = {10, 20, 50, 100, 200};
  int nN = 5;

  for (int ni = 0; ni < nN; ni++) {
    int N = N_vals[ni];
    double lam[MAXN], mu_arr[MAXN];
    for (int k = 0; k < N; k++)
      lam[k] = ZETA_ZEROS_2000[k];
    for (int k = 0; k < N - 1; k++)
      mu_arr[k] = 0.5 * (ZETA_ZEROS_2000[k] + ZETA_ZEROS_2000[k + 1]);

    double log_w[MAXN];
    int sign_w[MAXN];
    compute_raw_weights(lam, mu_arr, N, log_w, sign_w);

    double log_sum = logsumexp(log_w, N);
    double w_norm[MAXN];
    for (int k = 0; k < N; k++)
      w_norm[k] = exp(log_w[k] - log_sum);

    printf("  N = %d:\n", N);
    printf("  %4s  %10s  %5s  %12s  %12s  %12s  %12s\n", "j", "gamma_j", "eps",
           "Im(G_raw)*eps", "w_j_raw", "Im(G_norm)*eps", "w_j_norm");
    printf("  %4s  %10s  %5s  %12s  %12s  %12s  %12s\n", "---", "---", "---",
           "---", "---", "---", "---");

    for (int ji = 0; ji < nj; ji++) {
      int j = j_vals[ji];
      if (j >= N)
        continue;
      for (int ei = 0; ei < neps; ei++) {
        double eps = eps_vals[ei];
        double zr = ZETA_ZEROS_2000[j];

        Cpx gr = raw_cauchy_log(lam, log_w, sign_w, N, zr, eps);
        Cpx gn = norm_cauchy(lam, w_norm, N, zr, eps);

        double wj_raw = exp(log_w[j]);
        double wj_norm = w_norm[j];

        printf("  %4d  %10.4f  %5.3f  %12.4e  %12.4e  %12.4e  %12.4e\n", j, zr,
               eps, gr.im * eps, wj_raw, gn.im * eps, wj_norm);
      }
    }
    printf("\n");
  }
}

static void computation4(void) {
  printf("\n");
  printf("================================================================\n");
  printf("  COMPUTATION 4: KEY - does G_raw(z)/zeta'/zeta(z) converge?\n");
  printf(
      "================================================================\n\n");

  int mid_j[] = {4, 9, 14, 19, 29, 49};
  int nmid = 6;
  double zi = 0.1;

  int N_vals[] = {10, 20, 30, 50, 100, 200, 500, 1000};
  int nN = 8;

  double zr_pts[6];
  for (int mi = 0; mi < nmid; mi++)
    zr_pts[mi] =
        0.5 * (ZETA_ZEROS_2000[mid_j[mi]] + ZETA_ZEROS_2000[mid_j[mi] + 1]);

  printf("  Ratio R_N(z) = |G_raw(z)| / |zeta'/zeta(z)|\n\n");
  printf("  %6s", "N");
  for (int mi = 0; mi < nmid; mi++)
    printf("  z=%.1f", zr_pts[mi]);
  printf("  %12s  %12s\n", "mean(R)", "std(R)");
  printf("  %6s", "---");
  for (int mi = 0; mi < nmid; mi++)
    printf("  %12s", "---");
  printf("  %12s  %12s\n", "---", "---");

  for (int ni = 0; ni < nN; ni++) {
    int N = N_vals[ni];
    double lam[MAXN], mu_arr[MAXN];
    for (int k = 0; k < N; k++)
      lam[k] = ZETA_ZEROS_2000[k];
    for (int k = 0; k < N - 1; k++)
      mu_arr[k] = 0.5 * (ZETA_ZEROS_2000[k] + ZETA_ZEROS_2000[k + 1]);

    double log_w[MAXN];
    int sign_w[MAXN];
    compute_raw_weights(lam, mu_arr, N, log_w, sign_w);

    printf("  %6d", N);
    double sum_r = 0.0, sum_r2 = 0.0;
    int nr = 0;
    for (int mi = 0; mi < nmid; mi++) {
      if (mid_j[mi] >= N - 1) {
        printf("  %12s", "---");
        continue;
      }
      double zr = zr_pts[mi];
      Cpx gr = raw_cauchy_log(lam, log_w, sign_w, N, zr, zi);
      Cpx zpz = zeta_prime_over_zeta(zr, NPRIMES);
      double r = (cpx_abs(zpz) > 1e-15) ? cpx_abs(gr) / cpx_abs(zpz) : 0.0;
      printf("  %12.4e", r);
      sum_r += r;
      sum_r2 += r * r;
      nr++;
    }
    if (nr > 0) {
      double mean = sum_r / (double)nr;
      double var = sum_r2 / (double)nr - mean * mean;
      double std_val = (var > 0.0) ? sqrt(var) : 0.0;
      printf("  %12.4e  %12.4e", mean, std_val);
    }
    printf("\n");
  }
}

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

  printf("================================================================\n");
  printf(
      "  Unnormalized Cauchy Transform: G_raw(z) = sum w_k^raw/(gam_k - z)\n");
  printf("  vs zeta'/zeta(1/2 + iz)\n");
  printf("  Using %d zeta zeros, %d primes\n", N_REF_2000, NPRIMES);
  printf("  Weights computed in log-space for stability\n");
  printf("================================================================\n");

  computation1();
  computation2();
  computation3();
  computation4();

  printf("\nDone.\n");
  return 0;
}
