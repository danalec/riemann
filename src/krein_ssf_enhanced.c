/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Enhanced Krein SSF
 * @paper   yamaguchi-rh-2026.tex, §10.6
 * @theorem Theorem II
 * @proof   Obstruction analysis
 * @step    4
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
#define NPRIMES200 200

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
  double *w = (double *)xmalloc((size_t)N * sizeof(double));
  double ws = 0.0;
  for (int k = 0; k < N; k++) {
    double log_nm = 0.0;
    double log_dn = 0.0;
    int sign_nm = 1;
    int sign_dn = 1;
    for (int j = 0; j < N - 1; j++) {
      double d = lam[k] - mu[j];
      if (d == 0.0) {
        free(w);
        return -1;
      }
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
      if (d == 0.0) {
        free(w);
        return -1;
      }
      if (d < 0.0) {
        sign_dn = -sign_dn;
        d = -d;
      }
      log_dn += log(d);
    }
    double log_w = log_nm - log_dn;
    if (log_w > 700.0) {
      free(w);
      return -2;
    }
    w[k] = (double)(sign_nm * sign_dn) * exp(log_w);
    if (w[k] < 0.0) {
      free(w);
      return -2;
    }
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
  free(w);
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
    double im_sum = 0.0;
    Cpx m_cur = cpx_inv(cpx_make(a_arr[N - 1] - E, -eps));
    for (int j = N - 2; j >= 0; j--) {
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
  }
  return tr_im;
}

static Cpx zeta_prime_over_zeta_eps(double E, double eps_val, int np) {
  Cpx sum = cpx_make(0.0, 0.0);
  for (int pi = 0; pi < np; pi++) {
    double p = (double)PRIMES200[pi];
    double lp = log(p);
    double s_im = E + eps_val;
    Cpx phase = cpx_p_ik(p, s_im);
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

static void computation1(void) {
  printf("\n");
  printf("================================================================\n");
  printf("  COMPUTATION 1: Euler Product Convergence Study\n");
  printf(
      "================================================================\n\n");

  int zero_idx[] = {9, 19, 29, 49};
  int n_zero_idx = 4;

  printf("  1a: zeta'/zeta at zeta zeros z = gamma_k + i*epsilon\n");
  printf("  Using Euler form: p^{is} = e^{is ln p} via cpx_p_ik\n\n");

  double eps_vals[] = {0.1, 0.01, 0.001};
  int n_eps = 3;
  int prime_counts[] = {10, 20, 50, 100, 200};
  int n_pc = 5;

  for (int zi = 0; zi < n_zero_idx; zi++) {
    int k = zero_idx[zi];
    double gamma_k = ZETA_ZEROS_2000[k];
    printf("  --- At gamma_%d = %.6f ---\n", k + 1, gamma_k);
    printf("  %5s  %8s  %12s  %12s  %12s\n", "nprime", "eps", "Re(z'/z)",
           "Im(z'/z)", "|z'/z|");
    printf("  %5s  %8s  %12s  %12s  %12s\n", "---", "---", "---", "---", "---");

    for (int ei = 0; ei < n_eps; ei++) {
      double eps = eps_vals[ei];
      for (int pi = 0; pi < n_pc; pi++) {
        int np = prime_counts[pi];
        Cpx zpz = zeta_prime_over_zeta_eps(gamma_k, eps, np);
        printf("  %5d  %8.3f  %12.4f  %12.4f  %12.4f\n", np, eps, zpz.re,
               zpz.im, cpx_abs(zpz));
      }
      printf("\n");
    }
  }

  printf("  1b: zeta'/zeta at regular points (midpoints between zeros)\n\n");

  int mid_idx[] = {9, 19, 29, 49};
  int n_mid = 4;

  for (int mi = 0; mi < n_mid; mi++) {
    int k = mid_idx[mi];
    double E_mid = 0.5 * (ZETA_ZEROS_2000[k] + ZETA_ZEROS_2000[k + 1]);
    printf("  --- Midpoint gamma_%d..%d = %.6f ---\n", k + 1, k + 2, E_mid);
    printf("  %5s  %8s  %12s  %12s  %12s\n", "nprime", "eps", "Re(z'/z)",
           "Im(z'/z)", "|z'/z|");
    printf("  %5s  %8s  %12s  %12s  %12s\n", "---", "---", "---", "---", "---");

    for (int ei = 0; ei < n_eps; ei++) {
      double eps = eps_vals[ei];
      for (int pi = 0; pi < n_pc; pi++) {
        int np = prime_counts[pi];
        Cpx zpz = zeta_prime_over_zeta_eps(E_mid, eps, np);
        printf("  %5d  %8.3f  %12.6f  %12.6f  %12.6f\n", np, eps, zpz.re,
               zpz.im, cpx_abs(zpz));
      }
      printf("\n");
    }
  }
}

static void computation2(void) {
  printf("\n");
  printf("================================================================\n");
  printf("  COMPUTATION 2: Krein SSF Pointwise at Large N\n");
  printf(
      "================================================================\n\n");

  int N_vals[] = {50, 100, 200, 500};
  int n_N = 4;
  int np = NPRIMES200;

  for (int ni = 0; ni < n_N; ni++) {
    int N = N_vals[ni];
    printf("  --- N = %d ---\n", N);

    double *lam = (double *)xmalloc((size_t)N * sizeof(double));
    double *mu_arr = (double *)xmalloc((size_t)(N - 1) * sizeof(double));
    double *a_arr = (double *)xmalloc((size_t)N * sizeof(double));
    double *b_arr = (double *)xmalloc((size_t)(N - 1) * sizeof(double));

    for (int k = 0; k < N; k++)
      lam[k] = ZETA_ZEROS_2000[k];
    for (int k = 0; k < N - 1; k++)
      mu_arr[k] = 0.5 * (ZETA_ZEROS_2000[k] + ZETA_ZEROS_2000[k + 1]);

    if (deboor(lam, mu_arr, N, a_arr, b_arr) != 0) {
      printf("  dBG failed for N=%d\n\n", N);
      free(lam);
      free(mu_arr);
      free(a_arr);
      free(b_arr);
      continue;
    }

    printf("  xi(E) = SturmCount(E) - theta(E)/pi\n");
    printf("  S(E)  = -(1/pi) * arg zeta(1/2+iE)  [200 primes]\n\n");

    int nE = 40;
    double E_min = 10.0;
    double E_max = ZETA_ZEROS_2000[N - 1];
    double dE = (E_max - E_min) / (double)(nE - 1);

    printf("  %8s  %8s  %10s  %10s  %12s\n", "E", "N_zeta", "xi(E)", "S(E)",
           "|Delta|");
    printf("  %8s  %8s  %10s  %10s  %12s\n", "---", "---", "---", "---", "---");

    double sum_sq = 0.0;
    double max_delta = 0.0;
    int count = 0;

    for (int i = 0; i < nE; i++) {
      double E = E_min + (double)i * dE;
      int n_zeta = sturm_count(a_arr, b_arr, N, E);
      double n_weyl = theta_s(E) / M_PI;
      double xi = (double)n_zeta - n_weyl;

      double arg_z = arg_zeta_half(E, np);
      double S = arg_z / M_PI;

      double delta = fabs(xi - S);
      sum_sq += delta * delta;
      if (delta > max_delta)
        max_delta = delta;
      count++;

      printf("  %8.3f  %8d  %10.4f  %10.4f  %12.6f\n", E, n_zeta, xi, S, delta);
    }

    double rms = (count > 0) ? sqrt(sum_sq / (double)count) : 0.0;
    printf("\n  RMS |xi - S| = %.6f,  max |xi - S| = %.6f\n\n", rms, max_delta);

    free(lam);
    free(mu_arr);
    free(a_arr);
    free(b_arr);
  }
}

static void computation3(void) {
  printf("\n");
  printf("================================================================\n");
  printf("  COMPUTATION 3: Resolvent-Based SSF with Tighter epsilon\n");
  printf(
      "================================================================\n\n");

  int N_vals[] = {50, 100};
  int n_N = 2;
  double eps_vals[] = {0.1, 0.01, 0.001};
  int n_eps = 3;

  for (int ni = 0; ni < n_N; ni++) {
    int N = N_vals[ni];
    printf("  --- N = %d ---\n\n", N);

    double *lam = (double *)xmalloc((size_t)N * sizeof(double));
    double *mu_arr = (double *)xmalloc((size_t)(N - 1) * sizeof(double));
    double *a_arr = (double *)xmalloc((size_t)N * sizeof(double));
    double *b_arr = (double *)xmalloc((size_t)(N - 1) * sizeof(double));

    for (int k = 0; k < N; k++)
      lam[k] = ZETA_ZEROS_2000[k];
    for (int k = 0; k < N - 1; k++)
      mu_arr[k] = 0.5 * (ZETA_ZEROS_2000[k] + ZETA_ZEROS_2000[k + 1]);

    if (deboor(lam, mu_arr, N, a_arr, b_arr) != 0) {
      printf("  dBG failed for N=%d\n\n", N);
      free(lam);
      free(mu_arr);
      free(a_arr);
      free(b_arr);
      continue;
    }

    int nE = 30;
    double E_min = 10.0;
    double E_max = ZETA_ZEROS_2000[N - 1];
    double dE = (E_max - E_min) / (double)(nE - 1);

    for (int ei = 0; ei < n_eps; ei++) {
      double eps = eps_vals[ei];
      printf("  epsilon = %.3f:\n", eps);
      printf("  %8s  %10s  %12s  %12s\n", "E", "xi_sturm", "xi_res", "|Delta|");
      printf("  %8s  %10s  %12s  %12s\n", "---", "---", "---", "---");

      double sum_sq = 0.0;
      double max_d = 0.0;
      int cnt = 0;

      for (int i = 0; i < nE; i++) {
        double E = E_min + (double)i * dE;
        int n_zeta = sturm_count(a_arr, b_arr, N, E);
        double n_weyl = theta_s(E) / M_PI;
        double xi_sturm = (double)n_zeta - n_weyl;

        double tr_im = resolvent_trace_imag(a_arr, b_arr, N, E, eps);
        double xi_res = -tr_im / M_PI;

        double delta = fabs(xi_sturm - xi_res);
        sum_sq += delta * delta;
        if (delta > max_d)
          max_d = delta;
        cnt++;

        printf("  %8.3f  %10.4f  %12.4f  %12.6f\n", E, xi_sturm, xi_res, delta);
      }
      double rms = (cnt > 0) ? sqrt(sum_sq / (double)cnt) : 0.0;
      printf("  RMS=%.6f  max=%.6f\n\n", rms, max_d);
    }

    free(lam);
    free(mu_arr);
    free(a_arr);
    free(b_arr);
  }
}

static void computation4(void) {
  printf("\n");
  printf("================================================================\n");
  printf("  COMPUTATION 4: Spectral Measure Cumulative\n");
  printf(
      "================================================================\n\n");

  int N_vals[] = {100, 200};
  int n_N = 2;
  int np = NPRIMES200;

  for (int ni = 0; ni < n_N; ni++) {
    int N = N_vals[ni];
    printf("  --- N = %d ---\n\n", N);

    double *lam = (double *)xmalloc((size_t)N * sizeof(double));
    double *mu_arr = (double *)xmalloc((size_t)(N - 1) * sizeof(double));
    double *a_arr = (double *)xmalloc((size_t)N * sizeof(double));
    double *b_arr = (double *)xmalloc((size_t)(N - 1) * sizeof(double));

    for (int k = 0; k < N; k++)
      lam[k] = ZETA_ZEROS_2000[k];
    for (int k = 0; k < N - 1; k++)
      mu_arr[k] = 0.5 * (ZETA_ZEROS_2000[k] + ZETA_ZEROS_2000[k + 1]);

    if (deboor(lam, mu_arr, N, a_arr, b_arr) != 0) {
      printf("  dBG failed for N=%d\n\n", N);
      free(lam);
      free(mu_arr);
      free(a_arr);
      free(b_arr);
      continue;
    }

    double *weights = (double *)xmalloc((size_t)N * sizeof(double));
    {
      double *w = (double *)xmalloc((size_t)N * sizeof(double));
      double ws = 0.0;
      for (int k = 0; k < N; k++) {
        double nm = 1.0;
        for (int j = 0; j < N - 1; j++)
          nm *= lam[k] - mu_arr[j];
        double dn = 1.0;
        for (int j = 0; j < N; j++)
          if (j != k)
            dn *= lam[k] - lam[j];
        w[k] = nm / dn;
        if (w[k] < 0.0)
          w[k] = 0.0;
        ws += w[k];
      }
      for (int k = 0; k < N; k++)
        weights[k] = w[k] / ws;
      free(w);
    }

    printf(
        "  M(E) = Sum_{gamma_k <= E} w_k    (cumulative spectral measure)\n");
    printf("  N(E)/N = uniform counting\n");
    printf("  (1/pi)*arg zeta(1/2+iE) = zeta argument\n\n");

    int nE = 40;
    double E_min = ZETA_ZEROS_2000[0] - 2.0;
    double E_max = ZETA_ZEROS_2000[N - 1] + 2.0;
    double dE = (E_max - E_min) / (double)(nE - 1);

    printf("  %8s  %10s  %10s  %10s  %10s  %10s\n", "E", "M(E)", "N(E)/N",
           "arg/pi", "|M-NE|", "|M-arg|");
    printf("  %8s  %10s  %10s  %10s  %10s  %10s\n", "---", "---", "---", "---",
           "---", "---");

    double sum_sq_mn = 0.0, sum_sq_ma = 0.0;
    double max_mn = 0.0, max_ma = 0.0;
    int cnt = 0;

    for (int i = 0; i < nE; i++) {
      double E = E_min + (double)i * dE;

      double meas = 0.0;
      int count_below = 0;
      for (int k = 0; k < N; k++) {
        if (ZETA_ZEROS_2000[k] <= E) {
          meas += weights[k];
          count_below++;
        }
      }

      double NE_N = (double)count_below / (double)N;

      double arg_z = arg_zeta_half(E, np);
      double arg_pi = arg_z / M_PI;

      double delta_mn = fabs(meas - NE_N);
      double delta_ma = fabs(meas - arg_pi);
      sum_sq_mn += delta_mn * delta_mn;
      sum_sq_ma += delta_ma * delta_ma;
      if (delta_mn > max_mn)
        max_mn = delta_mn;
      if (delta_ma > max_ma)
        max_ma = delta_ma;
      cnt++;

      printf("  %8.3f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f\n", E, meas, NE_N,
             arg_pi, delta_mn, delta_ma);
    }

    double rms_mn = (cnt > 0) ? sqrt(sum_sq_mn / (double)cnt) : 0.0;
    double rms_ma = (cnt > 0) ? sqrt(sum_sq_ma / (double)cnt) : 0.0;
    printf("\n  |M(E) - N(E)/N|:  RMS = %.6f, max = %.6f\n", rms_mn, max_mn);
    printf("  |M(E) - arg/pi|:  RMS = %.6f, max = %.6f\n\n", rms_ma, max_ma);

    free(lam);
    free(mu_arr);
    free(a_arr);
    free(b_arr);
    free(weights);
  }
}

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

  printf("================================================================\n");
  printf("  Enhanced Krein Spectral Shift Function Comparison\n");
  printf("  2000 zeta zeros, 200 primes, tighter epsilon\n");
  printf("  Krein SSF identity: xi(E) = -(1/pi) arg zeta(1/2+iE)\n");
  printf("  All prime-frequency ops use Euler form via cpx_p_ik\n");
  printf("================================================================\n");

  computation1();
  computation2();
  computation3();
  computation4();

  printf("\nDone.\n");
  return 0;
}
