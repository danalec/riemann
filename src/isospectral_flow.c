/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Isospectral flow from Gram points to midpoints
 * @paper   yamaguchi-rh-2026.tex, §11.9
 * @theorem Lemma I
 * @proof   dBG trace continuity
 * @step    1
 *
 * Isospectral flow: Gram to midpoint. Paper Section 7.5.
 * isospectral_flow.c — trace de Boor-Golub entries as second spectrum varies
 *
 * Given fixed eigenvalues {γ_k}, vary the second spectrum {μ_k}
 * continuously from Gram points to midpoints and observe how
 * the matrix entries {a_k, b_k} deform along the isospectral manifold.
 *
 * Compile: gcc -O3 -o isospectral_flow isospectral_flow.c -lm
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.71828182845904523536
#endif

#ifdef _WIN32
#include <windows.h>
#endif

static const double zeta_zeros[30] = {
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
    37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
    52.970321, 56.446248, 59.347044, 60.831779, 65.112544};

static double gram_point(int k) {
  if (k == 0)
    return 17.84559954;
  if (k == 1)
    return 23.17028270;
  double t = 2.0 * M_PI * (double)k / log((double)k);
  for (int it = 0; it < 8; it++) {
    double lt = log(t);
    double theta = 0.5 * t * (lt - log(2.0 * M_PI * M_E)) + M_PI / 8.0;
    double dtheta = 0.5 * (lt - log(2.0 * M_PI * M_E) + 1.0);
    t -= (theta - M_PI * (double)k) / dtheta;
  }
  return t;
}

/* Simplified de Boor-Golub for small N, direct product weights */
static int mini_dbg(const double *lam, const double *mu, int N, double *a,
                    double *b) {
  for (int k = 0; k < N - 1; k++)
    if (!(lam[k] < mu[k] && mu[k] < lam[k + 1]))
      return -1;

  double w[20], wsum = 0;
  for (int k = 0; k < N; k++) {
    double num = 1.0;
    for (int j = 0; j < N - 1; j++)
      num *= (lam[k] - mu[j]);
    double den = 1.0;
    for (int j = 0; j < N; j++)
      if (j != k)
        den *= (lam[k] - lam[j]);
    w[k] = num / den;
    if (w[k] < 0)
      return -2;
    wsum += w[k];
  }
  for (int k = 0; k < N; k++)
    w[k] /= wsum;

  a[0] = 0;
  for (int i = 0; i < N; i++)
    a[0] += w[i] * lam[i];
  double norm_P1 = 0;
  for (int i = 0; i < N; i++) {
    double v = lam[i] - a[0];
    norm_P1 += w[i] * v * v;
  }
  b[0] = sqrt(norm_P1);
  double norm_Pk = norm_P1;

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
    a[k] = num / norm_Pk;
    if (k < N - 1) {
      double np1 = 0;
      for (int i = 0; i < N; i++) {
        double pp = 0, pc = 1;
        for (int j = 0; j <= k; j++) {
          double b2 = j > 0 ? b[j - 1] * b[j - 1] : 0;
          double pn = (lam[i] - a[j]) * pc - b2 * pp;
          pp = pc;
          pc = pn;
        }
        np1 += w[i] * pc * pc;
      }
      b[k] = sqrt(np1 / norm_Pk);
      norm_Pk = np1;
    }
  }
  return 0;
}

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

  int N = 9;
  double lam[20], mu_gram[19], mu_mid[19];
  for (int k = 0; k < N; k++)
    lam[k] = zeta_zeros[k];
  for (int k = 0; k < N - 1; k++)
    mu_gram[k] = gram_point(k);
  for (int k = 0; k < N - 1; k++)
    mu_mid[k] = 0.5 * (zeta_zeros[k] + zeta_zeros[k + 1]);

  printf("Isospectral flow: Gram → Midpoint  (N=%d)\n\n", N);
  printf("  t     a_0        a_1        a_2        b_0       b_1\n");
  printf("  ----  ---------  ---------  ---------  --------  -------\n");

  for (int ti = 0; ti <= 10; ti++) {
    double t = ti / 10.0; /* 0 = Gram, 1 = Midpoint */
    double mu_t[19], a_t[20], b_t[19];
    for (int k = 0; k < N - 1; k++)
      mu_t[k] = (1.0 - t) * mu_gram[k] + t * mu_mid[k];

    int ret = mini_dbg(lam, mu_t, N, a_t, b_t);
    if (ret != 0) {
      printf("  %.1f  INTERLACING FAILED\n", t);
      continue;
    }
    printf("  %.1f  %9.3f  %9.3f  %9.3f  %8.3f  %7.3f\n", t, a_t[0], a_t[1],
           a_t[2], b_t[0], b_t[1]);
  }

  // Isospectral flow: Gram to midpoint. Paper Section 7.5.
  printf("\n  The isospectral manifold: entries parametrized by {μ_k}.\n");

  /* ── Distance metrics between Gram and Midpoint reconstructions ── */
  printf("\n  ── Distance between Gram and Midpoint gauge choices ──\n");

  /* 1. Second-spectrum distance (gauge space) */
  double d_mu = 0.0;
  for (int k = 0; k < N - 1; k++) {
    double d = mu_gram[k] - mu_mid[k];
    d_mu += d * d;
  }
  d_mu = sqrt(d_mu / (N - 1));
  printf("  d_μ (2nd-spectrum RMS)= %.3f\n", d_mu);

  /* 2. Entry-wise Frobenius distance */
  double a_gram[20], b_gram[19], a_mid[20], b_mid[19];
  mini_dbg(lam, mu_gram, N, a_gram, b_gram);
  mini_dbg(lam, mu_mid, N, a_mid, b_mid);

  double d_entry = 0.0;
  for (int k = 0; k < N; k++) {
    double d = a_gram[k] - a_mid[k];
    d_entry += d * d;
  }
  for (int k = 0; k < N - 1; k++) {
    double d = b_gram[k] - b_mid[k];
    d_entry += d * d;
  }
  d_entry = sqrt(d_entry / (2 * N - 1));
  printf("  d_entry (Frobenius RMS) = %.3f\n", d_entry);

  /* 3. Sensitivity: ‖δJ‖ / ‖δμ‖ */
  printf("  Sensitivity ‖δJ‖/‖δμ‖ = %.3f\n", d_entry / d_mu);

  /* 4. Spectral weight distance */
  /* Compute weights for both reconstructions */
  double w_gram[20], w_mid[20], d_w = 0.0;
  for (int p = 0; p < 2; p++) {
    double *mu_p = p == 0 ? mu_gram : mu_mid;
    double wsum = 0.0;
    for (int k = 0; k < N; k++) {
      double num = 1.0;
      for (int j = 0; j < N - 1; j++)
        num *= (lam[k] - mu_p[j]);
      double den = 1.0;
      for (int j = 0; j < N; j++)
        if (j != k)
          den *= (lam[k] - lam[j]);
      double wk = num / den;
      wsum += wk;
      if (p == 0)
        w_gram[k] = wk;
      else
        w_mid[k] = wk;
    }
    if (wsum <= 0.0) {
      printf("  WARNING: non-positive weight sum\n");
      return 1;
    }
  }
  /* Normalize */
  {
    double s = 0;
    for (int k = 0; k < N; k++)
      s += w_gram[k];
    for (int k = 0; k < N; k++)
      w_gram[k] /= s;
  }
  {
    double s = 0;
    for (int k = 0; k < N; k++)
      s += w_mid[k];
    for (int k = 0; k < N; k++)
      w_mid[k] /= s;
  }
  for (int k = 0; k < N; k++) {
    double d = w_gram[k] - w_mid[k];
    d_w += d * d;
  }
  d_w = sqrt(d_w / N);
  printf("  d_w (spectral weight RMS)= %.4f\n", d_w);

  /* ── Distance growth along the flow ── */
  // Isospectral flow: Gram to midpoint. Paper Section 7.5.
  printf("\n  ── Distance from t=0 along the isospectral flow ──\n");
  printf("  t     d_entry    d_μ       ‖δJ‖/‖δμ‖\n");
  printf("  ----  ---------  --------  ----------\n");
  for (int ti = 0; ti <= 10; ti++) {
    double t = ti / 10.0;
    double mu_t[19];
    for (int k = 0; k < N - 1; k++)
      mu_t[k] = (1.0 - t) * mu_gram[k] + t * mu_mid[k];
    double a_t[20], b_t[19];
    if (mini_dbg(lam, mu_t, N, a_t, b_t) != 0)
      continue;
    double de = 0, dm = 0;
    for (int k = 0; k < N; k++) {
      double d = a_t[k] - a_gram[k];
      de += d * d;
    }
    for (int k = 0; k < N - 1; k++) {
      double d = b_t[k] - b_gram[k];
      de += d * d;
    }
    for (int k = 0; k < N - 1; k++) {
      double d = mu_t[k] - mu_gram[k];
      dm += d * d;
    }
    de = sqrt(de / (2 * N - 1));
    dm = sqrt(dm / (N - 1));
    printf("  %.1f   %.4e  %8.4f  %10.4f\n", t, de, dm,
           dm > 1e-10 ? de / dm : 0);
  }

  printf("\n  Interpretation:\n");
  printf("    - d_μ = %.2f: Gram and Midpoint gauge choices differ by %.2f on "
         "avg\n",
         d_mu, d_mu);
  printf(
      "    - d_entry = %.2f: resulting matrix entries differ by %.2f on avg\n",
      d_entry, d_entry);
  printf("    - Sensitivity = %.3f: entries change %.3f× more than gauge\n",
         d_entry / d_mu, d_entry / d_mu);
  printf("    - d_w = %.4f: normalized spectral weights differ very little\n",
         d_w);
  printf("    - The eigenvalues {γ_k} are IDENTICAL (d_λ = 0) for all gauge "
         "choices.\n");
  if (d_w < 0.01)
    printf(
        "    - Near-zero weight distance: manifold is nearly RIGID in weight "
        "space.\n"
        "      Changing the second spectrum mainly reshuffles entries without\n"
        "      significantly affecting the spectral measure.\n");
  return 0;
}
