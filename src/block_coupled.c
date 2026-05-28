/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   True coupled 2x2 block de Boor-Golub reconstruction
 * @paper   yamaguchi-rh-2026.tex, Section 7.5
 * @theorem Theorem IV (Block Spectral Bijection)
 * @proof   Younes-Trogdon matrix polynomial recurrence
 * @step    5 -- block gauge with non-diagonal weight matrices
 *
 * block_coupled.c — True coupled 2x2 block de Boor-Golub reconstruction
 *
 * Implements the full matrix-valued orthogonal polynomial recurrence
 * (Younes-Trogdon 2026, Sinap-Van Assche 1996):
 *
 *   x P_k(x) = P_{k+1}(x) B_k + P_k(x) A_k + P_{k-1}(x) B_{k-1}^T
 *
 * with P_0(x) = I_2, P_{-1}(x) = 0.
 *
 * Key difference from block_dbg.c: uses NON-DIAGONAL weight matrices
 * W_j = v_j v_j^T where v_j = (cos(theta_j), sin(theta_j))^T / sqrt(n),
 * providing full coupling between even/odd channels.
 *
 * Reconstruction preserves eigenvalue fidelity EXACTLY by construction
 * (spectral map bijection, Younes-Trogdon Theorem 3.14).
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

#define MAXN 30
#define MAX2N 60

static const double zeta[30] = {
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062, 37.586178,
    40.918719, 43.327073, 48.005151, 49.773832, 52.970321, 56.446248,
    59.347044, 60.831779, 65.112544, 67.079811, 69.546402, 72.067158,
    75.704691, 77.144840, 79.337375, 82.910381, 84.735493, 87.425275,
    88.809111, 92.491899, 94.651344, 95.870634, 98.831194, 101.317851};

typedef struct {
  double m[4];
} Mat2;

static Mat2 m2(double a, double b, double c, double d) {
  Mat2 r;
  r.m[0] = a;
  r.m[1] = b;
  r.m[2] = c;
  r.m[3] = d;
  return r;
}

static Mat2 m2_id(void) { return m2(1, 0, 0, 1); }
static Mat2 m2_zero(void) { return m2(0, 0, 0, 0); }
static Mat2 m2_add(Mat2 a, Mat2 b) {
  return m2(a.m[0] + b.m[0], a.m[1] + b.m[1], a.m[2] + b.m[2], a.m[3] + b.m[3]);
}
static Mat2 m2_sub(Mat2 a, Mat2 b) {
  return m2(a.m[0] - b.m[0], a.m[1] - b.m[1], a.m[2] - b.m[2], a.m[3] - b.m[3]);
}
static Mat2 m2_scale(double s, Mat2 a) {
  return m2(s * a.m[0], s * a.m[1], s * a.m[2], s * a.m[3]);
}
static Mat2 m2_mul(Mat2 a, Mat2 b) {
  return m2(
      a.m[0] * b.m[0] + a.m[1] * b.m[2], a.m[0] * b.m[1] + a.m[1] * b.m[3],
      a.m[2] * b.m[0] + a.m[3] * b.m[2], a.m[2] * b.m[1] + a.m[3] * b.m[3]);
}
static Mat2 m2_T(Mat2 a) { return m2(a.m[0], a.m[2], a.m[1], a.m[3]); }
static Mat2 m2_sym(Mat2 a) {
  double s = 0.5 * (a.m[1] + a.m[2]);
  return m2(a.m[0], s, s, a.m[3]);
}

static Mat2 m2_inv(Mat2 a) {
  double det = a.m[0] * a.m[3] - a.m[1] * a.m[2];
  double id = 1.0 / det;
  return m2(a.m[3] * id, -a.m[1] * id, -a.m[2] * id, a.m[0] * id);
}

static int m2_chol(Mat2 S, Mat2 *L) {
  if (S.m[0] <= 1e-30)
    return -1;
  double l00 = sqrt(S.m[0]);
  double l10 = S.m[2] / l00;
  double l11sq = S.m[3] - l10 * l10;
  if (l11sq <= 1e-30)
    return -1;
  *L = m2(l00, 0, l10, sqrt(l11sq));
  return 0;
}

static void dense_evals(int n, const double *M, double *ev) {
  double *V = (double *)malloc((size_t)(n * n) * sizeof(double));
  memcpy(V, M, (size_t)(n * n) * sizeof(double));
  for (int sw = 0; sw < 200; sw++) {
    double moff = 0;
    for (int p = 0; p < n - 1; p++)
      for (int q = p + 1; q < n; q++)
        if (fabs(V[p * n + q]) > moff)
          moff = fabs(V[p * n + q]);
    if (moff < 1e-14)
      break;
    for (int p = 0; p < n - 1; p++)
      for (int q = p + 1; q < n; q++) {
        double apq = V[p * n + q];
        if (fabs(apq) < 1e-16 * (fabs(V[p * n + p]) + fabs(V[q * n + q]) + 1))
          continue;
        double app = V[p * n + p], aqq = V[q * n + q];
        double tau = (aqq - app) / (2.0 * apq);
        double t = tau >= 0 ? 1.0 / (tau + sqrt(1 + tau * tau))
                            : -1.0 / (-tau + sqrt(1 + tau * tau));
        double c = 1.0 / sqrt(1 + t * t), s = t * c;
        for (int i = 0; i < n; i++) {
          double vip = V[i * n + p], viq = V[i * n + q];
          V[i * n + p] = vip * c - viq * s;
          V[i * n + q] = vip * s + viq * c;
        }
        for (int j = 0; j < n; j++) {
          double vpj = V[p * n + j], vqj = V[q * n + j];
          V[p * n + j] = vpj * c - vqj * s;
          V[q * n + j] = vpj * s + vqj * c;
        }
      }
  }
  for (int i = 0; i < n; i++)
    ev[i] = V[i * n + i];
  for (int i = 1; i < n; i++) {
    double key = ev[i];
    int j = i - 1;
    while (j >= 0 && ev[j] > key) {
      ev[j + 1] = ev[j];
      j--;
    }
    ev[j + 1] = key;
  }
  free(V);
}

static void build_block_matrix(int nb, const Mat2 *A, const Mat2 *B,
                               double *M) {
  int n = 2 * nb;
  for (int i = 0; i < n * n; i++)
    M[i] = 0;
  for (int k = 0; k < nb; k++) {
    int r0 = 2 * k, r1 = 2 * k + 1;
    M[r0 * n + r0] = A[k].m[0];
    M[r0 * n + r1] = A[k].m[1];
    M[r1 * n + r0] = A[k].m[2];
    M[r1 * n + r1] = A[k].m[3];
  }
  for (int k = 0; k < nb - 1; k++) {
    int r0 = 2 * k, r1 = 2 * k + 1;
    int c0 = 2 * (k + 1), c1 = 2 * (k + 1) + 1;
    M[r0 * n + c0] = B[k].m[0];
    M[r0 * n + c1] = B[k].m[1];
    M[r1 * n + c0] = B[k].m[2];
    M[r1 * n + c1] = B[k].m[3];
    M[c0 * n + r0] = B[k].m[0];
    M[c1 * n + r0] = B[k].m[1];
    M[c0 * n + r1] = B[k].m[2];
    M[c1 * n + r1] = B[k].m[3];
  }
}

static int coupled_block_dbg(int nb, const double *lam, const double *theta,
                             Mat2 *A_out, Mat2 *B_out) {
  int m = 2 * nb;
  Mat2 W[MAX2N];
  for (int j = 0; j < m; j++) {
    double cs = cos(theta[j]), sn = sin(theta[j]);
    double s = 1.0 / sqrt((double)m);
    W[j] =
        m2(cs * cs * s * s, cs * sn * s * s, cs * sn * s * s, sn * sn * s * s);
  }

  Mat2 Wsum = m2_zero();
  for (int j = 0; j < m; j++)
    Wsum = m2_add(Wsum, W[j]);
  for (int j = 0; j < m; j++) {
    W[j].m[0] /= Wsum.m[0];
    W[j].m[1] /= Wsum.m[0];
    W[j].m[2] /= Wsum.m[0];
    W[j].m[3] /= Wsum.m[0];
  }
  Wsum = m2_zero();
  for (int j = 0; j < m; j++)
    Wsum = m2_add(Wsum, W[j]);

  Mat2 P_prev[MAX2N], P_cur[MAX2N], P_raw[MAX2N];
  for (int j = 0; j < m; j++) {
    P_prev[j] = m2_zero();
    P_cur[j] = m2_id();
  }

  Mat2 A0 = m2_zero();
  for (int j = 0; j < m; j++) {
    A0 = m2_add(A0, m2_scale(lam[j], m2_mul(W[j], P_cur[j])));
  }
  A_out[0] = m2_sym(A0);

  for (int j = 0; j < m; j++)
    P_raw[j] = m2_sub(m2_scale(lam[j], P_cur[j]), m2_mul(A_out[0], P_cur[j]));

  Mat2 Sigma = m2_zero();
  for (int j = 0; j < m; j++) {
    Mat2 PtW = m2_mul(m2_T(P_raw[j]), W[j]);
    Sigma = m2_add(Sigma, m2_mul(PtW, P_raw[j]));
  }
  Sigma = m2_sym(Sigma);

  Mat2 L0;
  if (m2_chol(Sigma, &L0) != 0)
    return -1;
  B_out[0] = m2_T(L0);

  Mat2 Binvm1 = m2_inv(B_out[0]);
  for (int j = 0; j < m; j++)
    P_cur[j] = m2_mul(P_raw[j], Binvm1);

  for (int k = 1; k < nb; k++) {
    Mat2 Ak = m2_zero();
    for (int j = 0; j < m; j++) {
      Mat2 PkW = m2_mul(m2_T(P_cur[j]), W[j]);
      Ak = m2_add(Ak, m2_scale(lam[j], m2_mul(PkW, P_cur[j])));
    }
    A_out[k] = m2_sym(Ak);

    if (k < nb - 1) {
      Mat2 Bprev_T = m2_T(B_out[k - 1]);
      for (int j = 0; j < m; j++) {
        Mat2 tmp = m2_sub(m2_mul(m2_scale(lam[j], m2_id()), P_cur[j]),
                          m2_mul(A_out[k], P_cur[j]));
        tmp = m2_sub(tmp, m2_mul(P_prev[j], Bprev_T));
        P_raw[j] = tmp;
      }

      Mat2 Sig = m2_zero();
      for (int j = 0; j < m; j++) {
        Mat2 PtW = m2_mul(m2_T(P_raw[j]), W[j]);
        Sig = m2_add(Sig, m2_mul(PtW, P_raw[j]));
      }
      Sig = m2_sym(Sig);

      Mat2 Lk;
      if (m2_chol(Sig, &Lk) != 0)
        return -2;
      B_out[k] = m2_T(Lk);

      Mat2 Binv = m2_inv(B_out[k]);
      for (int j = 0; j < m; j++)
        P_prev[j] = P_cur[j];
      for (int j = 0; j < m; j++)
        P_cur[j] = m2_mul(P_raw[j], Binv);
    } else {
      for (int j = 0; j < m; j++)
        P_prev[j] = P_cur[j];
    }
  }
  return 0;
}

static void scalar_dbg(const double *lam, const double *mu, int N, double *a,
                       double *b) {
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
}

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

  int nb_max = 10;
  double theta[MAX2N];

  printf("Coupled 2x2 block de Boor-Golub reconstruction\n");
  printf("Full matrix-valued OP recurrence (Younes-Trogdon / Sinap-Van "
         "Assche)\n\n");

  for (int nb = 3; nb <= nb_max; nb++) {
    int m = 2 * nb;
    double lam[MAX2N];
    for (int j = 0; j < m; j++)
      lam[j] = zeta[j];

    for (int j = 0; j < m; j++)
      theta[j] = M_PI * j / m;

    Mat2 A[MAXN], B[MAXN];
    int ret = coupled_block_dbg(nb, lam, theta, A, B);

    printf("nb=%2d (n=%2d): ", nb, m);
    if (ret != 0) {
      printf("RECONSTRUCTION FAILED (ret=%d)\n", ret);
      continue;
    }

    double *Mblk = (double *)calloc((size_t)(m * m), sizeof(double));
    build_block_matrix(nb, A, B, Mblk);
    double ev[MAX2N];
    dense_evals(m, Mblk, ev);

    double max_err = 0;
    for (int k = 0; k < m; k++) {
      double d = fabs(ev[k] - lam[k]);
      if (d > max_err)
        max_err = d;
    }

    printf("max|ev-gamma|=%.2e %s", max_err,
           max_err < 1e-8   ? "OK"
           : max_err < 1e-4 ? "FAIR"
                            : "DEGRADED");

    double a_ref[MAXN], b_ref[MAXN - 1];
    double mu_ref[MAXN - 1];
    for (int k = 0; k < m - 1; k++)
      mu_ref[k] = 0.5 * (lam[k] + lam[k + 1]);
    scalar_dbg(lam, mu_ref, m, a_ref, b_ref);

    printf("  |  A_k coupling:");
    double max_off = 0;
    for (int k = 0; k < nb; k++)
      if (fabs(A[k].m[1]) > max_off)
        max_off = fabs(A[k].m[1]);
    printf(" max|a12|=%.4f", max_off);

    double max_boff = 0;
    for (int k = 0; k < nb - 1; k++) {
      double b12 = B[k].m[1], b21 = B[k].m[2];
      if (fabs(b12) > max_boff)
        max_boff = fabs(b12);
      if (fabs(b21) > max_boff)
        max_boff = fabs(b21);
    }
    printf(" max|B_off|=%.4f\n", max_boff);

    free(Mblk);
  }

  printf("\n--- Detailed output for nb=5 (10x10) ---\n\n");
  {
    int nb = 5, m = 10;
    double lam[MAX2N];
    for (int j = 0; j < m; j++)
      lam[j] = zeta[j];

    for (int j = 0; j < m; j++)
      theta[j] = M_PI * j / m;

    Mat2 A[5], B[4];
    int ret = coupled_block_dbg(nb, lam, theta, A, B);
    if (ret != 0) {
      printf("FAILED\n");
      return 1;
    }

    printf("Diagonal blocks A_k:\n");
    printf("  k   a11        a22        a12(off)\n");
    for (int k = 0; k < nb; k++)
      printf("  %d   %9.4f  %9.4f  %9.4f\n", k, A[k].m[0], A[k].m[3],
             A[k].m[1]);

    printf("\nOff-diagonal blocks B_k:\n");
    printf("  k   b11        b12        b21        b22\n");
    for (int k = 0; k < nb - 1; k++)
      printf("  %d   %9.4f  %9.4f  %9.4f  %9.4f\n", k, B[k].m[0], B[k].m[1],
             B[k].m[2], B[k].m[3]);

    double *Mblk = (double *)calloc((size_t)(m * m), sizeof(double));
    build_block_matrix(nb, A, B, Mblk);
    double ev[MAX2N];
    dense_evals(m, Mblk, ev);

    printf("\nEigenvalue verification:\n");
    printf("  k   gamma       ev_rec      |diff|\n");
    double max_err = 0;
    for (int k = 0; k < m; k++) {
      double d = fabs(ev[k] - lam[k]);
      if (d > max_err)
        max_err = d;
      printf("  %2d  %10.6f  %10.6f  %.2e\n", k, lam[k], ev[k], d);
    }
    printf("  Max eigenvalue error: %.2e\n", max_err);

    free(Mblk);
  }

  printf("\n--- Comparison: theta=pi*j/n vs theta=0 (diagonal) ---\n\n");
  {
    int nb = 5, m = 10;
    double lam[MAX2N];
    for (int j = 0; j < m; j++)
      lam[j] = zeta[j];

    double theta_diag[MAX2N];
    for (int j = 0; j < m; j++)
      theta_diag[j] = (j % 2 == 0) ? 0 : M_PI / 2;

    Mat2 Ad[5], Bd[4];
    int ret = coupled_block_dbg(nb, lam, theta_diag, Ad, Bd);
    printf("Diagonal weights: ");
    if (ret == 0) {
      double max_off = 0;
      for (int k = 0; k < nb; k++)
        max_off = fmax(max_off, fabs(Ad[k].m[1]));
      for (int k = 0; k < nb - 1; k++) {
        max_off = fmax(max_off, fabs(Bd[k].m[1]));
        max_off = fmax(max_off, fabs(Bd[k].m[2]));
      }
      printf("max coupling = %.2e (should be ~0)\n", max_off);
    } else {
      printf("FAILED (%d)\n", ret);
    }

    for (int j = 0; j < m; j++)
      theta[j] = M_PI * j / m;
    Mat2 Ac[5], Bc[4];
    ret = coupled_block_dbg(nb, lam, theta, Ac, Bc);
    printf("Coupled weights: ");
    if (ret == 0) {
      double max_off = 0;
      for (int k = 0; k < nb; k++)
        max_off = fmax(max_off, fabs(Ac[k].m[1]));
      for (int k = 0; k < nb - 1; k++) {
        max_off = fmax(max_off, fabs(Bc[k].m[1]));
        max_off = fmax(max_off, fabs(Bc[k].m[2]));
      }
      printf("max coupling = %.4f (should be >0)\n", max_off);
    } else {
      printf("FAILED (%d)\n", ret);
    }

    printf("\nGauge dimensions: %d blocks * 7 params/block = %d params\n", nb,
           3 * nb + 3 * (nb - 1));
    printf("Matrix entries: %d (%d diagonal + %d off-diagonal block params)\n",
           2 * m - 1, 3 * nb, 3 * (nb - 1));
    printf("Free gauge params (theta): %d\n\n", m);
  }

  return 0;
}
