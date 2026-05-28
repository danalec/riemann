/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Optimize isospectral gauge to match prime perturbation
 * @paper   yamaguchi-rh-2026.tex, Section 7.3
 * @theorem Lemma III (Isospectral Gauge Freedom)
 * @proof   Gradient descent on isospectral manifold for prime correlation
 * @step    4 -- gauge optimization for prime pattern matching
 *
 * prime_fit.c — optimize isospectral gauge to match prime perturbation
 *
 * Given eigenvalues {γ_k}, search the isospectral manifold (parametrized
 * by {μ_k}) to find the Jacobi matrix whose entries best match the
 * prime-correlated perturbation pattern:
 *   δa_k = Σ_p α_p · sin(log p · k)
 *   δb_k = Σ_p α_p · cos(log p · k)    with α_p = -log(p)/(2π√p)
 *
 * Compile: gcc -O3 -o prime_fit prime_fit.c -lm
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

static const double zeta[30] = {
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062, 37.586178, 40.918719,
    43.327073, 48.005151, 49.773832, 52.970321, 56.446248, 59.347044, 60.831779,
    65.112544, 67.079811, 69.546402, 72.067158, 75.704691, 77.144840, 79.337375,
    82.910381, 84.735493, 87.425275, 88.809111};

static int deboor(const double *lam, const double *mu, int N, double *a,
                  double *b) {
  for (int k = 0; k < N - 1; k++)
    if (!(lam[k] < mu[k] && mu[k] < lam[k + 1]))
      return -1;
  double w[50], ws = 0;
  for (int k = 0; k < N; k++) {
    double n = 1;
    for (int j = 0; j < N - 1; j++)
      n *= lam[k] - mu[j];
    double d = 1;
    for (int j = 0; j < N; j++)
      if (j != k)
        d *= lam[k] - lam[j];
    w[k] = n / d;
    if (w[k] < 0)
      return -2;
    ws += w[k];
  }
  for (int k = 0; k < N; k++)
    w[k] /= ws;
  a[0] = 0;
  for (int i = 0; i < N; i++)
    a[0] += w[i] * lam[i];
  double np1 = 0;
  for (int i = 0; i < N; i++) {
    double v = lam[i] - a[0];
    np1 += w[i] * v * v;
  }
  b[0] = sqrt(np1);
  double npk = np1;
  for (int k = 1; k < N; k++) {
    double num = 0;
    for (int i = 0; i < N; i++) {
      double pp = 0, pc = 1;
      for (int j = 0; j < k; j++) {
        double b2 = j > 0 ? b[j - 1] * b[j - 1] : 0;
        double pn = (lam[i] - a[j]) * pc - b2 * pp;
        pp = pc;
        pc = pn;
      }
      num += w[i] * lam[i] * pc * pc;
    }
    a[k] = num / npk;
    if (k < N - 1) {
      double npx = 0;
      for (int i = 0; i < N; i++) {
        double pp = 0, pc = 1;
        for (int j = 0; j <= k; j++) {
          double b2 = j > 0 ? b[j - 1] * b[j - 1] : 0;
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

/* Target prime perturbation */
static void prime_target(int N, double *da, double *db) {
  int primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
  for (int k = 0; k < N; k++) {
    da[k] = 0;
    if (k < N - 1)
      db[k] = 0;
    for (int pi = 0; pi < 10; pi++) {
      double w = log(primes[pi]);
      double a = -log(primes[pi]) / (2 * M_PI * sqrt(primes[pi]));
      da[k] += a * sin(w * k);
      if (k < N - 1)
        db[k] += a * cos(w * k);
    }
  }
}

/* Distance between reconstructed entries and reference + prime target */
static double fit_error(const double *a, const double *b, const double *aref,
                        const double *bref, int N, double scale) {
  double da_target[50], db_target[49];
  prime_target(N, da_target, db_target);
  double err = 0;
  for (int k = 0; k < N; k++) {
    double expected = aref[k] + scale * da_target[k];
    double d = a[k] - expected;
    err += d * d;
  }
  for (int k = 0; k < N - 1; k++) {
    double expected = bref[k] + scale * db_target[k];
    double d = b[k] - expected;
    err += d * d;
  }
  return sqrt(err / (2 * N - 1));
}

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif
  int N = 9;
  double lam[30], mu_mid[29];
  for (int k = 0; k < N; k++)
    lam[k] = zeta[k];
  for (int k = 0; k < N - 1; k++)
    mu_mid[k] = 0.5 * (zeta[k] + zeta[k + 1]);

  /* Reference: free Jacobi reconstructed from midpoint spectrum */
  double a_ref[30], b_ref[29];
  deboor(lam, mu_mid, N, a_ref, b_ref);

  /* Grid search: vary μ_0 and μ_1, find best fit to prime pattern */
  double mu0_mid = mu_mid[0], mu1_mid = mu_mid[1];
  double best_err = 1e300, best_mu0 = mu0_mid, best_mu1 = mu1_mid;

  printf("Prime fit optimization on isospectral manifold (N=%d)\n\n", N);
  printf("  mu_0       mu_1       fit_err    a_0        b_0\n");
  printf("  ---------  ---------  ---------  ---------  ---------\n");

  for (int i = 0; i <= 20; i++) {
    double mu0 = mu0_mid + (i - 10) * 0.1;
    if (!(lam[0] < mu0 && mu0 < lam[1]))
      continue;
    for (int j = 0; j <= 20; j++) {
      double mu1 = mu1_mid + (j - 10) * 0.1;
      if (!(lam[1] < mu1 && mu1 < lam[2]))
        continue;

      double mu[29];
      memcpy(mu, mu_mid, (size_t)(N - 1) * sizeof(double));
      mu[0] = mu0;
      mu[1] = mu1;

      double a[30], b[29];
      if (deboor(lam, mu, N, a, b) != 0)
        continue;
      double err = fit_error(a, b, a_ref, b_ref, N, 1.0);
      printf("  %9.3f  %9.3f  %9.3f  %9.3f  %9.3f\n", mu0, mu1, err, a[0],
             b[0]);
      if (err < best_err) {
        best_err = err;
        best_mu0 = mu0;
        best_mu1 = mu1;
      }
    }
  }

  printf("\n  Best gauge (2-param): mu_0=%.3f mu_1=%.3f (error=%.3f)\n",
         best_mu0, best_mu1, best_err);
  printf("  Midpoint gauge: mu_0=%.3f mu_1=%.3f\n", mu0_mid, mu1_mid);

  /* ── Full 8-parameter Nelder-Mead optimization ── */
  printf("\n  ═══════════════════════════════════════════\n");
  printf("  Full Nelder-Mead: optimizing all %d mu parameters\n", N - 1);
  printf("  ═══════════════════════════════════════════\n\n");

  int dim = N - 1;
  double (*simplex)[29] =
      (double (*)[29])malloc((size_t)(dim + 1) * sizeof(double[29]));
  double *fv = (double *)malloc((size_t)(dim + 1) * sizeof(double));
  double *centroid = (double *)malloc((size_t)dim * sizeof(double));
  double *xr = (double *)malloc((size_t)dim * sizeof(double));
  double *xe = (double *)malloc((size_t)dim * sizeof(double));
  double *xc = (double *)malloc((size_t)dim * sizeof(double));

  /* Initialize simplex around midpoint gauge */
  for (int i = 0; i <= dim; i++) {
    for (int d = 0; d < dim; d++) {
      simplex[i][d] = mu_mid[d];
      if (i > 0 && i - 1 == d)
        simplex[i][d] += 0.5; /* spread simplex */
    }
    /* Clamp to interlacing bounds */
    for (int d = 0; d < dim; d++) {
      if (simplex[i][d] <= lam[d])
        simplex[i][d] = lam[d] + 0.01;
      if (simplex[i][d] >= lam[d + 1])
        simplex[i][d] = lam[d + 1] - 0.01;
    }
    /* Evaluate */
    double a_t[30], b_t[29];
    if (deboor(lam, simplex[i], N, a_t, b_t) != 0)
      fv[i] = 1e300;
    else
      fv[i] = fit_error(a_t, b_t, a_ref, b_ref, N, 1.0);
  }

  printf("  iter  best_err\n");
  printf("  ----  ---------\n");
  int max_iter = 200;
  for (int iter = 0; iter < max_iter; iter++) {
    /* Find hi, lo, second-hi */
    int hi = 0, lo = 0;
    for (int i = 1; i <= dim; i++) {
      if (fv[i] > fv[hi])
        hi = i;
      if (fv[i] < fv[lo])
        lo = i;
    }
    int s_hi = (hi == 0) ? 1 : 0;
    for (int i = 0; i <= dim; i++) {
      if (i == hi)
        continue;
      if (fv[i] > fv[s_hi])
        s_hi = i;
    }

    if (iter % 20 == 0)
      printf("  %4d  %.4f\n", iter, fv[lo]);

    /* Centroid (excluding hi) */
    for (int d = 0; d < dim; d++) {
      double sum = 0;
      for (int i = 0; i <= dim; i++)
        if (i != hi)
          sum += simplex[i][d];
      centroid[d] = sum / (double)dim;
    }

    /* Reflection */
    int valid = 1;
    for (int d = 0; d < dim; d++) {
      xr[d] = 2 * centroid[d] - simplex[hi][d];
      if (!(lam[d] < xr[d] && xr[d] < lam[d + 1]))
        valid = 0;
    }
    double fr = 1e300;
    if (valid) {
      double a_t[30], b_t[29];
      if (deboor(lam, xr, N, a_t, b_t) == 0)
        fr = fit_error(a_t, b_t, a_ref, b_ref, N, 1.0);
    }

    if (fr < fv[lo]) {
      /* Expansion */
      for (int d = 0; d < dim; d++) {
        xe[d] = 3 * centroid[d] - 2 * simplex[hi][d];
        if (!(lam[d] < xe[d] && xe[d] < lam[d + 1])) {
          xe[d] = xr[d];
        }
      }
      double fe = 1e300;
      valid = 1;
      for (int d = 0; d < dim; d++)
        if (!(lam[d] < xe[d] && xe[d] < lam[d + 1]))
          valid = 0;
      if (valid) {
        double a_t[30], b_t[29];
        if (deboor(lam, xe, N, a_t, b_t) == 0)
          fe = fit_error(a_t, b_t, a_ref, b_ref, N, 1.0);
      }
      if (fe < fr) {
        for (int d = 0; d < dim; d++)
          simplex[hi][d] = xe[d];
        fv[hi] = fe;
      } else {
        for (int d = 0; d < dim; d++)
          simplex[hi][d] = xr[d];
        fv[hi] = fr;
      }
    } else if (fr < fv[s_hi]) {
      for (int d = 0; d < dim; d++)
        simplex[hi][d] = xr[d];
      fv[hi] = fr;
    } else {
      /* Contract */
      if (fr < fv[hi]) {
        for (int d = 0; d < dim; d++)
          simplex[hi][d] = xr[d];
        fv[hi] = fr;
      }
      for (int d = 0; d < dim; d++) {
        xc[d] = 0.5 * (simplex[hi][d] + centroid[d]);
        if (!(lam[d] < xc[d] && xc[d] < lam[d + 1]))
          xc[d] = centroid[d];
      }
      double fc = 1e300;
      valid = 1;
      for (int d = 0; d < dim; d++)
        if (!(lam[d] < xc[d] && xc[d] < lam[d + 1]))
          valid = 0;
      if (valid) {
        double a_t[30], b_t[29];
        if (deboor(lam, xc, N, a_t, b_t) == 0)
          fc = fit_error(a_t, b_t, a_ref, b_ref, N, 1.0);
      }
      if (fc < fv[hi]) {
        for (int d = 0; d < dim; d++)
          simplex[hi][d] = xc[d];
        fv[hi] = fc;
      } else {
        /* Shrink */
        for (int i = 0; i <= dim; i++) {
          if (i == lo)
            continue;
          for (int d = 0; d < dim; d++) {
            simplex[i][d] = 0.5 * (simplex[i][d] + simplex[lo][d]);
            if (simplex[i][d] <= lam[d])
              simplex[i][d] = lam[d] + 0.01;
            if (simplex[i][d] >= lam[d + 1])
              simplex[i][d] = lam[d + 1] - 0.01;
          }
          double a_t[30], b_t[29];
          if (deboor(lam, simplex[i], N, a_t, b_t) != 0)
            fv[i] = 1e300;
          else
            fv[i] = fit_error(a_t, b_t, a_ref, b_ref, N, 1.0);
        }
      }
    }
  }

  /* Best result */
  int best_idx = 0;
  for (int i = 1; i <= dim; i++)
    if (fv[i] < fv[best_idx])
      best_idx = i;
  printf("  %4d  %.4f (final)\n", max_iter, fv[best_idx]);

  printf("\n  ── Best NM gauge (all 8 params) vs Midpoint ──\n");
  printf("  k   mu_best    mu_mid     δμ\n");
  printf("  ---  ---------  ---------  -------\n");
  for (int k = 0; k < dim; k++)
    printf("  %2d   %9.3f  %9.3f  %+7.3f\n", k, simplex[best_idx][k], mu_mid[k],
           simplex[best_idx][k] - mu_mid[k]);

  /* Compare entries */
  double a_best[30], b_best[29];
  deboor(lam, simplex[best_idx], N, a_best, b_best);
  printf("\n  ── Best entries vs Midpoint entries ──\n");
  printf("  k   a_best     a_mid      δa        b_best     b_mid      δb\n");
  printf("  ---  ---------  ---------  --------  ---------  ---------  "
         "--------\n");
  for (int k = 0; k < N; k++) {
    double da = a_best[k] - a_ref[k];
    double db = k < N - 1 ? b_best[k] - b_ref[k] : 0;
    printf("  %2d   %9.3f  %9.3f  %+8.3f  %9.3f  %9.3f  %+8.3f\n", k, a_best[k],
           a_ref[k], da, k < N - 1 ? b_best[k] : 0.0,
           k < N - 1 ? b_ref[k] : 0.0, db);
  }

  /* Compare to prime target */
  printf("\n  ── NM best entries vs Midpoint (norm) ──\n");
  double err_nm = fit_error(a_best, b_best, a_ref, b_ref, N, 1.0);
  double a_mid[30], b_mid[29];
  deboor(lam, mu_mid, N, a_mid, b_mid);
  double err_mid = fit_error(a_mid, b_mid, a_ref, b_ref, N, 1.0);
  printf("  NM best error: %.4f\n", err_nm);
  printf("  Midpoint error: %.4f\n", err_mid);
  printf("  Improvement: %.1f%%\n", 100.0 * (err_mid - err_nm) / err_mid);

  free(simplex);
  free(fv);
  free(centroid);
  free(xr);
  free(xe);
  free(xc);
  return 0;
}
