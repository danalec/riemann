/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   m-zeta bridge: dispersion relation
 * @paper   yamaguchi-rh-2026.tex, §10.6
 * @theorem Theorem III
 * @proof   Level 3 barrier
 * @step    4
 *
 * m_zeta_bridge.c -- m(z) = zeta'/zeta(1/2+iz) bridge test
 *
 * THE BRIDGE IDENTITY (Level 3 barrier):
 *   m(z) = zeta'/zeta(1/2 + iz)    for all z in C \ R
 *
 * where m(z) = <e_0, (J - z)^{-1} e_0> is the Weyl-Titchmarsh m-function
 * of the fully corrected Gram Jacobi matrix.
 *
 * APPROACH: Dispersion relation on m/m_free
 *
 * From L8: arg m_0^free(E+i0) = -theta(E) + C_0 + O(1/E)  [verified]
 * For full: arg m_0(E+i0)     = -theta(E) - pi*S(E) + C   [from SSF]
 *
 * Therefore: arg(m_0/m_0^free) = -pi*S(E) + (C - C_0)     on real axis
 *
 * By dispersion relation (Herglotz representation):
 *   log(m_0/m_0^free)(z) = (1/pi) integral Im[log(m/m_free)](E)/(E-z) dE
 *                          + polynomial terms
 *
 * If Im[log(m/m_free)](E) = -pi*S(E) + const, then:
 *   (m_0/m_free)(z) = K * exp(-i*pi*S(z))   [analytic continuation]
 *
 * We test this by computing m_0(z) at COMPLEX z (away from poles)
 * and comparing to zeta'/zeta(1/2+iz).
 *
 * Compile: gcc -Wall -Wextra -Werror -O3 -Isrc
 *          -o m_zeta_bridge src/m_zeta_bridge.c -lm
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

/* ---- Complex arithmetic ---- */

typedef struct {
  double re, im;
} Cplx;

static Cplx cplx(double r, double i) {
  Cplx z = {r, i};
  return z;
}
static Cplx cadd(Cplx a, Cplx b) { return cplx(a.re + b.re, a.im + b.im); }
static Cplx csub(Cplx a, Cplx b) { return cplx(a.re - b.re, a.im - b.im); }
static Cplx cmul(Cplx a, Cplx b) {
  return cplx(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}
static Cplx cdiv(Cplx a, Cplx b) {
  double d = b.re * b.re + b.im * b.im;
  return cplx((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d);
}
static double cmod(Cplx z) { return sqrt(z.re * z.re + z.im * z.im); }
static Cplx zlog(Cplx z) { return cplx(log(cmod(z)), atan2(z.im, z.re)); }
static Cplx zexp(Cplx z) {
  double e = exp(z.re);
  return cplx(e * cos(z.im), e * sin(z.im));
}

/* ---- log-Gamma via Stirling with reflection ---- */

static Cplx lgamma_c(Cplx z);

static Cplx lgamma_c(Cplx z) {
  if (z.re < 0.5) {
    Cplx pz = cplx(M_PI * z.re, M_PI * z.im);
    Cplx spz = cplx(sin(pz.re) * cosh(pz.im), cos(pz.re) * sinh(pz.im));
    Cplx omz = cplx(1.0 - z.re, -z.im);
    return csub(cplx(log(M_PI), 0.0), cadd(zlog(spz), lgamma_c(omz)));
  }
  Cplx lz = zlog(z);
  Cplx zmh = csub(z, cplx(0.5, 0.0));
  Cplx r = csub(cmul(zmh, lz), z);
  r.re += 0.5 * log(2.0 * M_PI);
  double z2 = z.re * z.re + z.im * z.im;
  r = cadd(r, cplx(z.re / (12.0 * z2), -z.im / (12.0 * z2)));
  double z4 = z2 * z2;
  double z3re = z.re * z.re * z.re - 3.0 * z.re * z.im * z.im;
  double z3im = 3.0 * z.re * z.re * z.im - z.im * z.im * z.im;
  r = cadd(r, cplx(-z3re / (360.0 * z4 * z2), -z3im / (360.0 * z4 * z2)));
  return r;
}

/* ---- Zeta via approximate functional equation ---- */

static Cplx zeta_func(Cplx s) {
  double sigma = s.re, t = s.im;
  double at = fabs(t);

  if (at < 2.0) {
    if (sigma > 1.5) {
      double sum = 0.0;
      for (int n = 1; n <= 100000; n++)
        sum += pow((double)n, -sigma);
      return cplx(sum, 0.0);
    }
    double sum = 0.0;
    for (int n = 1; n <= 100000; n++) {
      double ns = pow((double)n, -sigma);
      sum += ns;
    }
    double cor = pow(100000.0, 1.0 - sigma) / (sigma - 1.0);
    return cplx(sum + cor, 0.0);
  }

  int M = (int)(sqrt(at / (2.0 * M_PI)) + 0.5);
  if (M < 1)
    M = 1;
  if (M > 5000)
    M = 5000;

  Cplx s1 = cplx(0.0, 0.0);
  for (int n = 1; n <= M; n++) {
    double ln = log((double)n);
    double f = pow((double)n, -sigma);
    s1 = cadd(s1, cplx(f * cos(-t * ln), f * sin(-t * ln)));
  }

  Cplx s2 = cplx((1.0 - s.re) * 0.5, -s.im * 0.5);
  Cplx sh = cplx(s.re * 0.5, s.im * 0.5);
  Cplx chi_log = csub(csub(lgamma_c(s2), lgamma_c(sh)),
                      cplx((s.re - 0.5) * log(M_PI), s.im * log(M_PI)));
  Cplx chi = zexp(chi_log);

  Cplx s2b = cplx(0.0, 0.0);
  double sm1 = sigma - 1.0;
  for (int n = 1; n <= M; n++) {
    double ln = log((double)n);
    double f = pow((double)n, sm1);
    s2b = cadd(s2b, cplx(f * cos(t * ln), f * sin(t * ln)));
  }

  return cadd(s1, cmul(chi, s2b));
}

/* ---- zeta'/zeta via numerical differentiation ---- */

static Cplx zeta_log_deriv(Cplx s) {
  double h = 1e-8;
  Cplx sp = cplx(s.re + h, s.im);
  Cplx z1 = zeta_func(s);
  Cplx z2 = zeta_func(sp);
  double mod1 = cmod(z1);
  double mod2 = cmod(z2);
  if (mod1 < 1e-30 || mod2 < 1e-30)
    return cplx(0.0, 0.0);
  Cplx dz = cdiv(csub(z2, z1), cplx(h, 0.0));
  return cdiv(dz, z1);
}

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
    /* Correction from R-vM */
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

/* ---- m-function via continued fraction ----
 *
 * m(z) = -1/(z - a_0 - b_0^2/(z - a_1 - b_1^2/(z - a_2 - ...)))
 *
 * Computed bottom-up: start with m_N = 0, recurse upward.
 * For large N, m_N -> 0 is a good approximation for Im(z) > 0.
 */

static Cplx m_function(int N, const double *a, const double *b, double zr,
                       double zi) {
  /* Start from bottom */
  Cplx m = cplx(0.0, 0.0);

  for (int k = N - 1; k >= 0; k--) {
    /* m_k = -1 / (z - a_k - b_k^2 / m_{k+1}) */
    double zma_re = zr - a[k];
    double zma_im = zi;
    if (k < N - 1) {
      double b2 = b[k] * b[k];
      double mag2 = m.re * m.re + m.im * m.im;
      if (mag2 > 1e-30) {
        double inv_re = m.re / mag2;
        double inv_im = -m.im / mag2;
        zma_re -= b2 * inv_re;
        zma_im -= b2 * inv_im;
      }
    }
    double mag2 = zma_re * zma_re + zma_im * zma_im;
    if (mag2 > 1e-30) {
      m = cplx(-zma_re / mag2, zma_im / mag2);
    }
  }
  return m;
}

/* ---- Free m-function (uncorrected Gram Jacobi) ---- */

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

/* ---- S(E) = N(E) - theta(E)/pi - 1 ---- */

static double S_of_E(double E) {
  /* At zeta zeros: S(gamma_k) = k - 1 - theta(gamma_k)/pi */
  /* Interpolate from reference data */
  if (E <= ZETA_ZEROS[0])
    return 0.0;
  int k = 0;
  while (k < N_REF - 1 && ZETA_ZEROS[k + 1] < E)
    k++;
  if (k >= N_REF - 1) {
    /* Extrapolate: S(T) ~ (1/pi) * (T log(T/2pi) - T)/... */
    return 0.0;
  }
  /* Linear interpolation */
  double gk = ZETA_ZEROS[k], gkp1 = ZETA_ZEROS[k + 1];
  double Sk = (double)k - theta_s(gk) / M_PI;
  double Skp1 = (double)(k + 1) - theta_s(gkp1) / M_PI;
  double t = (E - gk) / (gkp1 - gk);
  return (1.0 - t) * Sk + t * Skp1;
}

/* ---- Main tests ---- */

static void test_m_vs_zeta_direct(void) {
  printf("================================================================\n");
  printf("  TEST 1: m(z) vs zeta'/zeta(1/2+iz) -- Direct Comparison\n");
  printf("  At complex z (Im z > 0), away from poles\n");
  printf(
      "================================================================\n\n");

  int N = 100;
  double *a = (double *)malloc((size_t)N * sizeof(double));
  double *b = (double *)malloc((size_t)N * sizeof(double));
  double *af = (double *)malloc((size_t)N * sizeof(double));
  double *bf = (double *)malloc((size_t)N * sizeof(double));

  build_jacobi(N, a, b);
  build_free_jacobi(N, af, bf);

  /* Test points: z = zr + i*zi with zi > 0 (upper half-plane) */
  double zr_vals[] = {10.0, 30.0, 50.0, 80.0, 100.0};
  double zi_vals[] = {2.0, 5.0, 10.0, 20.0, 30.0};
  int nzr = 5, nzi = 5;

  printf("  %8s  %10s  %14s  %14s  %10s\n", "z", "s=1/2+iz", "|m_J|", "|z'/z|",
         "ratio");
  printf("  %8s  %10s  %14s  %14s  %10s\n", "------", "------", "------",
         "------", "------");

  for (int ir = 0; ir < nzr; ir++) {
    for (int ii = 0; ii < nzi; ii++) {
      double zr = zr_vals[ir];
      double zi = zi_vals[ii];

      /* s = 1/2 + iz = (1/2 - zi) + i*zr */
      double sr = 0.5 - zi;
      double si = zr;

      Cplx mj = m_function(N, a, b, zr, zi);
      Cplx mf = m_function(N, af, bf, zr, zi);
      Cplx s = cplx(sr, si);
      Cplx zld = zeta_log_deriv(s);

      double mj_mod = cmod(mj);
      (void)cmod(mf); /* unused */
      double zld_mod = cmod(zld);

      /* The bridge identity: m_J(z) = zeta'/zeta(1/2+iz) */
      /* But m_J has a different normalization: m_J ~ -1/z for large z */
      /* while zeta'/zeta ~ 0.5 log z. So we compare the RATIO */
      double ratio = (zld_mod > 1e-15) ? mj_mod / zld_mod : 0.0;

      printf("  %5.0f%+4.0fi  %5.0f%+5.0fi  %14.6e  %14.6e  %10.4f\n", zr, zi,
             sr, si, mj_mod, zld_mod, ratio);
    }
  }

  printf("\n  Note: m_J and zeta'/zeta have different asymptotics:\n");
  printf("  m_J(z) ~ -1/z  (free m-function)\n");
  printf("  zeta'/zeta(s) ~ (1/2) log(s/2pi)  (Stirling)\n");
  printf("  The BRIDGE is: m_J/m_free = (zeta'/zeta) / (m_free_analytic)\n");
  printf("  i.e., the RATIO of full to free should match.\n\n");

  /* Compare ratios: m_J/m_free vs (zeta'/zeta) / (analytic m_free) */
  printf("  Ratio comparison: (m_J/m_free) vs predicted\n\n");
  printf("  %8s  %12s  %12s  %10s  %10s  %10s\n", "z", "|m_J/m_f|", "S(E) est",
         "-pi*S", "arg(m_J/m_f)", "match?");
  printf("  %8s  %12s  %12s  %10s  %10s  %10s\n", "------", "--------",
         "--------", "------", "----------", "------");

  for (int ir = 0; ir < nzr; ir++) {
    double zr = zr_vals[ir];
    double zi = 1.0; /* Small imaginary part for boundary value */

    Cplx mj = m_function(N, a, b, zr, zi);
    Cplx mf = m_function(N, af, bf, zr, zi);
    Cplx ratio_c = cdiv(mj, mf);

    double S_E = S_of_E(zr);
    double neg_pi_S = -M_PI * S_E;
    double arg_ratio = atan2(ratio_c.im, ratio_c.re);

    /* Wrap arg to [-pi, pi] */
    while (arg_ratio > M_PI)
      arg_ratio -= 2.0 * M_PI;
    while (arg_ratio < -M_PI)
      arg_ratio += 2.0 * M_PI;

    double diff = fabs(arg_ratio - neg_pi_S);
    while (diff > M_PI)
      diff = 2.0 * M_PI - diff;

    printf("  %8.1f  %12.6f  %12.6f  %10.4f  %10.4f  %s\n", zr, cmod(ratio_c),
           S_E, neg_pi_S, arg_ratio, diff < 0.5 ? "CLOSE" : "FAR");
  }

  free(a);
  free(b);
  free(af);
  free(bf);
}

static void test_dispersion_relation(void) {
  printf(
      "\n================================================================\n");
  printf("  TEST 2: Dispersion Relation Test\n");
  printf("  log(m_J/m_free)(z) = (1/pi) int Im[log(m/m_free)](E)/(E-z) dE\n");
  printf("  If Im[log(m/m_free)](E) = -pi*S(E), then:\n");
  printf("  log(m_J/m_free)(z) = -int S(E)/(E-z) dE\n");
  printf(
      "================================================================\n\n");

  int N = 100;
  double *a = (double *)malloc((size_t)N * sizeof(double));
  double *b = (double *)malloc((size_t)N * sizeof(double));
  double *af = (double *)malloc((size_t)N * sizeof(double));
  double *bf = (double *)malloc((size_t)N * sizeof(double));

  build_jacobi(N, a, b);
  build_free_jacobi(N, af, bf);

  double z_test[][2] = {
      {50.0, 5.0},
      {50.0, 10.0},
      {100.0, 5.0},
      {100.0, 10.0},
  };
  int nzp = 4;

  for (int iz = 0; iz < nzp; iz++) {
    double zr = z_test[iz][0];
    double zi = z_test[iz][1];

    printf("  z = %.1f + %.1fi\n\n", zr, zi);

    /* Direct: log(m_J/m_free)(z) */
    Cplx mj = m_function(N, a, b, zr, zi);
    Cplx mf = m_function(N, af, bf, zr, zi);
    Cplx ratio = cdiv(mj, mf);
    Cplx log_ratio = zlog(ratio);

    printf("  Direct: log(m_J/m_free)(z) = %.6f + %.6fi\n", log_ratio.re,
           log_ratio.im);

    /* Dispersion: (1/pi) int S(E)/(E-z) dE */
    /* Numerical integration from E_min to E_max */
    double E_min = 10.0, E_max = 200.0;
    int npts = 10000;
    double dE = (E_max - E_min) / (double)(npts - 1);

    double disp_re = 0.0, disp_im = 0.0;
    for (int i = 0; i < npts; i++) {
      double E = E_min + i * dE;
      double S_E = S_of_E(E);
      double w = (i == 0 || i == npts - 1) ? 0.5 : 1.0;

      /* 1/(E - z) = (E - zr - i*zi) / ((E-zr)^2 + zi^2) */
      double den = (E - zr) * (E - zr) + zi * zi;
      double kernel_re = (E - zr) / den;
      double kernel_im = zi / den;

      disp_re += w * S_E * kernel_re;
      disp_im += w * S_E * kernel_im;
    }
    disp_re *= dE / M_PI;
    disp_im *= dE / M_PI;

    printf("  Dispersion: -int S/(E-z) dE = %.6f + %.6fi\n", -disp_re,
           -disp_im);

    printf("  Diff: Re = %.6f, Im = %.6f\n\n", fabs(log_ratio.re + disp_re),
           fabs(log_ratio.im + disp_im));
  }

  free(a);
  free(b);
  free(af);
  free(bf);
}

static void test_convergence_with_N(void) {
  printf(
      "\n================================================================\n");
  printf("  TEST 3: Convergence of m_J/m_free with N\n");
  printf("  At fixed z, does m_J/m_free stabilize as N increases?\n");
  printf(
      "================================================================\n\n");

  double zr = 50.0, zi = 5.0;
  int Ns[] = {25, 50, 100, 150, 200};
  int nN = 5;

  printf("  z = %.1f + %.1fi\n\n", zr, zi);
  printf("  %4s  %14s  %14s  %14s  %14s  %12s\n", "N", "m_J", "m_free",
         "m_J/m_free", "|m_J/m_free|", "arg(m_J/m_f)");
  printf("  %4s  %14s  %14s  %14s  %14s  %12s\n", "---", "---", "---", "---",
         "---", "---");

  for (int iN = 0; iN < nN; iN++) {
    int N = Ns[iN];
    double *a = (double *)malloc((size_t)N * sizeof(double));
    double *b = (double *)malloc((size_t)N * sizeof(double));
    double *af = (double *)malloc((size_t)N * sizeof(double));
    double *bf = (double *)malloc((size_t)N * sizeof(double));

    build_jacobi(N, a, b);
    build_free_jacobi(N, af, bf);

    Cplx mj = m_function(N, a, b, zr, zi);
    Cplx mf = m_function(N, af, bf, zr, zi);
    Cplx ratio = cdiv(mj, mf);
    double ratio_mod = cmod(ratio);
    double ratio_arg = atan2(ratio.im, ratio.re);

    printf("  %4d  %7.4f%+7.4fi  %7.4f%+7.4fi  %7.4f%+7.4fi  %14.6f  %12.4f\n",
           N, mj.re, mj.im, mf.re, mf.im, ratio.re, ratio.im, ratio_mod,
           ratio_arg);

    free(a);
    free(b);
    free(af);
    free(bf);
  }
}

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

  printf("################################################################\n");
  printf("#  m(z) = zeta'/zeta(1/2+iz) BRIDGE TEST\n");
  printf("#\n");
  printf("#  The Level 3 barrier: proving the m-function identity\n");
  printf("#  which is equivalent to the Riemann Hypothesis.\n");
  printf("#\n");
  printf("#  Strategy: Dispersion relation on m_J/m_free\n");
  printf("#  If Im[log(m_J/m_free)](E) = -pi*S(E) on real axis,\n");
  printf("#  then log(m_J/m_free)(z) = -int S(E)/(E-z) dE by dispersion.\n");
  printf("#  This would give m_J = c*m_free*exp(-i*pi*S(z)).\n");
  printf(
      "################################################################\n\n");

  test_m_vs_zeta_direct();
  test_dispersion_relation();
  test_convergence_with_N();

  printf("################################################################\n");
  printf("#  INTERPRETATION\n");
  printf("#\n");
  printf("#  Test 1: If m_J and zeta'/zeta have similar moduli and\n");
  printf("#  arg(m_J/m_free) approx -pi*S(E), the bridge is supported.\n");
  printf("#\n");
  printf("#  Test 2: If the dispersion integral matches log(m_J/m_free)\n");
  printf("#  at complex z, the Herglotz representation holds.\n");
  printf("#\n");
  printf("#  Test 3: If m_J/m_free stabilizes with N, the limit exists.\n");
  printf("################################################################\n");

  return 0;
}
