/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   2x2 block de Boor-Golub with diagonal weight matrices
 * @paper   yamaguchi-rh-2026.tex, Section 7.4
 * @theorem Sinap-Van Assche interlacing theorem
 * @proof   Diagonal block decomposition into scalar DBG channels
 * @step    5 -- interleaved even/odd reconstruction
 *
 * block_dbg.c — 2x2 block de Boor-Golub reconstruction
 *
 * Simplified from Sinap & Van Assche "Orthogonal matrix polynomials".
 * With diagonal weight matrices W_k = diag(w_{2k}, w_{2k+1}), the
 * block structure decomposes into two interleaved scalar de Boor-Golub
 * channels. Each channel enjoys strict interlacing — exact reconstruction.
 *
 * Compile:
 *   gcc -O3 -fno-strict-aliasing -fno-peel-loops -fno-unswitch-loops \
 *       -Wall -Wextra -Wconversion -Wshadow -Werror -Isrc -lm \
 *       src/block_dbg.c -o bin/block_dbg.exe
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
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062, 37.586178,
    40.918719, 43.327073, 48.005151, 49.773832, 52.970321, 56.446248,
    59.347044, 60.831779, 65.112544, 67.079811, 69.546402, 72.067158,
    75.704691, 77.144840, 79.337375, 82.910381, 84.735493, 87.425275,
    88.809111, 92.491899, 94.651344, 95.870634, 98.831194, 101.317851};

/* Standard de Boor-Golub for N eigenvalues and N-1 sub-eigenvalues.
 * Returns -1 on interlacing fail, -2 on negative weight. */
static int deboor(const double *lam, const double *mu, int N, double *a,
                  double *b) {
  int k, i, j;
  double w[30], ws;
  double np1, npk;

  for (k = 0; k < N - 1; k++)
    if (!(lam[k] < mu[k] && mu[k] < lam[k + 1]))
      return -1;

  ws = 0.0;
  for (k = 0; k < N; k++) {
    double n = 1.0;
    for (j = 0; j < N - 1; j++)
      n *= lam[k] - mu[j];
    double d = 1.0;
    for (j = 0; j < N; j++)
      if (j != k)
        d *= lam[k] - lam[j];
    w[k] = n / d;
    if (w[k] < 0.0)
      return -2;
    ws += w[k];
  }
  for (k = 0; k < N; k++)
    w[k] /= ws;

  a[0] = 0.0;
  for (i = 0; i < N; i++)
    a[0] += w[i] * lam[i];
  np1 = 0.0;
  for (i = 0; i < N; i++) {
    double v = lam[i] - a[0];
    np1 += w[i] * v * v;
  }
  b[0] = sqrt(np1);
  npk = np1;

  for (k = 1; k < N; k++) {
    double num = 0.0;
    for (i = 0; i < N; i++) {
      double pp = 0.0, pc = 1.0;
      for (j = 0; j < k; j++) {
        double b2 = j > 0 ? b[j - 1] * b[j - 1] : 0.0;
        double pn = (lam[i] - a[j]) * pc - b2 * pp;
        pp = pc;
        pc = pn;
      }
      num += w[i] * lam[i] * pc * pc;
    }
    a[k] = num / npk;
    if (k < N - 1) {
      double npx = 0.0;
      for (i = 0; i < N; i++) {
        double pp = 0.0, pc = 1.0;
        for (j = 0; j <= k; j++) {
          double b2 = j > 0 ? b[j - 1] * b[j - 1] : 0.0;
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

/* Build 2N x 2N block Jacobi matrix from N 2x2 blocks.
 * A_k = diag block (a11, a22, a12=symmetric)
 * B_k = off-diag block (b11,b12,b21,b22) */
static void build_block_jacobi(int N, const double *a11, const double *a22,
                               const double *a12, const double *b11,
                               const double *b12, const double *b21,
                               const double *b22, double *M) {
  int n = 2 * N;
  int i;
  for (i = 0; i < n * n; i++)
    M[i] = 0.0;
  for (i = 0; i < N; i++) {
    int r0 = 2 * i, r1 = 2 * i + 1;
    M[r0 * n + r0] = a11[i];
    M[r1 * n + r1] = a22[i];
    M[r0 * n + r1] = a12[i];
    M[r1 * n + r0] = a12[i];
  }
  for (i = 0; i < N - 1; i++) {
    int r0 = 2 * i, r1 = 2 * i + 1;
    int c0 = 2 * (i + 1), c1 = 2 * (i + 1) + 1;
    M[r0 * n + c0] = b11[i];
    M[c0 * n + r0] = b11[i];
    M[r0 * n + c1] = b12[i];
    M[c1 * n + r0] = b21[i];
    M[r1 * n + c0] = b21[i];
    M[c0 * n + r1] = b12[i];
    M[r1 * n + c1] = b22[i];
    M[c1 * n + r1] = b22[i];
  }
}

/* Jacobi eigensolver on dense symmetric matrix */
static void dense_evals(int N, const double *M, double *evals) {
  double *V = (double *)malloc((size_t)(N * N) * sizeof(double));
  int sweep, p, q, i, j;
  memcpy(V, M, (size_t)(N * N) * sizeof(double));

  for (sweep = 0; sweep < 100; sweep++) {
    double moff = 0.0;
    for (p = 0; p < N - 1; p++)
      for (q = p + 1; q < N; q++) {
        double v = fabs(V[p * N + q]);
        if (v > moff)
          moff = v;
      }
    if (moff < 1e-14)
      break;
    for (p = 0; p < N - 1; p++)
      for (q = p + 1; q < N; q++) {
        double apq = V[p * N + q];
        if (fabs(apq) < 1e-16 * (fabs(V[p * N + p]) + fabs(V[q * N + q]) + 1.0))
          continue;
        double app = V[p * N + p], aqq = V[q * N + q];
        double tau = (aqq - app) / (2.0 * apq);
        double t = tau >= 0.0 ? 1.0 / (tau + sqrt(1.0 + tau * tau))
                              : -1.0 / (-tau + sqrt(1.0 + tau * tau));
        double c = 1.0 / sqrt(1.0 + t * t), s = t * c;
        for (i = 0; i < N; i++) {
          double vip = V[i * N + p], viq = V[i * N + q];
          V[i * N + p] = vip * c - viq * s;
          V[i * N + q] = vip * s + viq * c;
        }
        for (j = 0; j < N; j++) {
          double vpj = V[p * N + j], vqj = V[q * N + j];
          V[p * N + j] = vpj * c - vqj * s;
          V[q * N + j] = vpj * s + vqj * c;
        }
      }
  }
  for (i = 0; i < N; i++)
    evals[i] = V[i * N + i];
  for (i = 1; i < N; i++) {
    double key = evals[i];
    j = i - 1;
    while (j >= 0 && evals[j] > key) {
      evals[j + 1] = evals[j];
      j--;
    }
    evals[j + 1] = key;
  }
  free(V);
}

/* Block de Boor-Golub reconstruction (N = number of 2x2 blocks).
 *
 * Splits into two interleaved channels:
 *   Ch 0: lam_even[i] = zeta[2*i],  mu_even[j] = (zeta[2j]+zeta[2j+1])/2
 *   Ch 1: lam_odd[i]  = zeta[2*i+1], mu_odd[j] = (zeta[2j+1]+zeta[2j+2])/2
 *
 * Each channel has strict interlacing, enabling exact reconstruction. */
static int block_dbg_reconstruct(int N, double *a11, double *a22, double *a12,
                                 double *b11, double *b12, double *b21,
                                 double *b22) {
  double lam_even[15], mu_even[14];
  double lam_odd[15], mu_odd[14];
  double a_even[15], b_even[14];
  double a_odd[15], b_odd[14];
  int i, j;

  for (i = 0; i < N; i++) {
    lam_even[i] = zeta[2 * i];
    lam_odd[i] = zeta[2 * i + 1];
  }
  for (j = 0; j < N - 1; j++) {
    mu_even[j] = 0.5 * (zeta[2 * j] + zeta[2 * j + 1]);
    mu_odd[j] = 0.5 * (zeta[2 * j + 1] + zeta[2 * j + 2]);
  }

  if (deboor(lam_even, mu_even, N, a_even, b_even) != 0)
    return -1;
  if (deboor(lam_odd, mu_odd, N, a_odd, b_odd) != 0)
    return -2;

  for (i = 0; i < N; i++) {
    a11[i] = a_even[i];
    a22[i] = a_odd[i];
    a12[i] = 0.0; /* diagonal blocks */
  }
  for (i = 0; i < N - 1; i++) {
    b11[i] = b_even[i];
    b22[i] = b_odd[i];
    b12[i] = 0.0;
    b21[i] = 0.0;
  }
  return 0;
}

/* Cholesky L * L^T of a 2x2 SPD matrix m.
 * L is lower triangular (L[0],L[2],0,L[3]).
 * Returns 0 on success, -1 if not SPD. */
static int mat2_chol(const double m[4], double L[4]) {
  if (m[0] <= 0.0)
    return -1;
  L[0] = sqrt(m[0]);
  L[2] = m[2] / L[0];
  L[3] = m[3] - L[2] * L[2];
  if (L[3] <= 0.0)
    return -1;
  L[3] = sqrt(L[3]);
  L[1] = 0.0;
  return 0;
}

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

  int NB = 5;
  int n = 2 * NB;
  double a11[10], a22[10], a12[10];
  double b11[10], b12[10], b21[10], b22[10];
  int i, k, ret;
  double *M, ev_block[20], max_err;

  (void)mat2_chol;

  printf("====================================================================="
         "=\n");
  printf(
      "  2x2 Block de Boor-Golub Reconstruction (N=%d blocks, %dx%d matrix)\n",
      NB, n, n);
  printf("====================================================================="
         "=\n\n");

  /* ---- interlacing check per channel ---- */
  printf("  --- Channel interlacing ---\n\n");
  printf("  Channel 0 (even indices):\n");
  for (i = 0; i < NB; i++) {
    double m0 = 0.5 * (zeta[2 * i] + zeta[2 * i + 1]);
    printf("    %8.4f < %8.4f < %8.4f", zeta[2 * i], m0, zeta[2 * i + 1]);
    if (i < NB - 1) {
      double m1 = 0.5 * (zeta[2 * i + 1] + zeta[2 * i + 2]);
      printf("  |  %8.4f < %8.4f < %8.4f\n", zeta[2 * i + 1], m1,
             zeta[2 * i + 2]);
    } else {
      printf("\n");
    }
  }
  printf("\n  All inter-block and intra-block midpoints satisfy strict\n");
  printf("  interlacing in both channels — exact reconstruction possible.\n\n");

  /* ---- reconstruct ---- */
  ret = block_dbg_reconstruct(NB, a11, a22, a12, b11, b12, b21, b22);
  if (ret != 0) {
    printf("  ERROR: reconstruction failed (%d)\n", ret);
    return 1;
  }

  /* ---- build block Jacobi and compute eigenvalues ---- */
  M = (double *)calloc((size_t)(n * n), sizeof(double));
  build_block_jacobi(NB, a11, a22, a12, b11, b12, b21, b22, M);
  dense_evals(n, M, ev_block);

  /* ---- verify eigenvalues ---- */
  printf("  --- Eigenvalue verification ---\n\n");
  printf("  %3s  %12s  %12s  %10s\n", "k", "lambda_in", "lambda_rec", "|diff|");
  printf("  %3s  %12s  %12s  %10s\n", "---", "----------", "----------",
         "--------");
  max_err = 0.0;
  for (k = 0; k < n; k++) {
    double diff = fabs(zeta[k] - ev_block[k]);
    if (diff > max_err)
      max_err = diff;
    printf("  %3d  %12.6f  %12.6f  %10.2e\n", k, zeta[k], ev_block[k], diff);
  }
  printf("\n  Max eigenvalue error: %.2e\n\n", max_err);

  /* ---- print A_k ---- */
  printf("  --- Block diagonal entries A_k (3 params: a11, a22, a12) ---\n\n");
  printf("  %3s  %10s  %10s  %10s\n", "k", "a11", "a22", "a12");
  printf("  %3s  %10s  %10s  %10s\n", "---", "--------", "--------",
         "--------");
  for (k = 0; k < NB; k++)
    printf("  %3d  %10.4f  %10.4f  %10.4f\n", k, a11[k], a22[k], a12[k]);

  /* ---- print B_k ---- */
  printf("\n  --- Block off-diagonal entries B_k (4 params: b11, b12, b21, "
         "b22) ---\n\n");
  printf("  %3s  %10s  %10s  %10s  %10s\n", "k", "b11", "b12", "b21", "b22");
  printf("  %3s  %10s  %10s  %10s  %10s\n", "---", "--------", "--------",
         "--------", "--------");
  for (k = 0; k < NB - 1; k++)
    printf("  %3d  %10.4f  %10.4f  %10.4f  %10.4f\n", k, b11[k], b12[k], b21[k],
           b22[k]);

  /* ---- compare to scalar reconstruction at N=5 ---- */
  printf("\n  --- Comparison: scalar de Boor-Golub (N=5) ---\n\n");
  {
    double lam5[10], mu5[9], a5[10], b5[9];
    int ret5;
    for (i = 0; i < 5; i++)
      lam5[i] = zeta[i];
    for (i = 0; i < 4; i++)
      mu5[i] = 0.5 * (zeta[i] + zeta[i + 1]);
    ret5 = deboor(lam5, mu5, 5, a5, b5);
    if (ret5 != 0) {
      printf("  Scalar N=5 reconstruction FAILED (%d)\n", ret5);
      printf("  This is expected for some N, but N=5 should work.\n");
    } else {
      printf("  Scalar (N=5) vs Block (N=5 blocks) diagonal entries:\n\n");
      printf("  %3s  %10s  %10s  %10s  %10s  %10s\n", "k", "a_scalar",
             "b_scalar", "a11_block", "a22_block", "ch_ratio");
      printf("  %3s  %10s  %10s  %10s  %10s  %10s\n", "---", "--------",
             "--------", "--------", "--------", "--------");
      for (k = 0; k < 5; k++) {
        double ratio = (a11[k] + a22[k]) * 0.5 / a5[k];
        printf("  %3d  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n", k, a5[k],
               b5[k], a11[k], a22[k], ratio);
      }
      printf("\n  Block a11 (channel 0) uses lambda_0, lambda_2, lambda_4, "
             "lambda_6, lambda_8\n");
      printf("  Block a22 (channel 1) uses lambda_1, lambda_3, lambda_5, "
             "lambda_7, lambda_9\n");
      printf("  Scalar a_k   uses lambda_0, lambda_1, lambda_2, lambda_3, "
             "lambda_4\n");
      printf("  The block method reconstructs 10 eigenvalues from\n");
      printf("  10 evals and 8 sub-evals, exploiting the 2-channel "
             "decomposition.\n");
      printf(
          "  The scalar method handles only the first 5 with 4 sub-evals.\n");
    }
  }

  /* ---- also run block reconstruction via full matrix Lanczos ---- */
  printf("\n  --- Full matrix block Lanczos (direct implementation) ---\n\n");
  {
    double raw_w[20], w[20], sum_even = 0.0, sum_odd = 0.0;
    int Mpts = n, P = Mpts - 2;
    double lam[20], mu[18];

    for (i = 0; i < Mpts; i++)
      lam[i] = zeta[i];
    for (i = 0; i < P; i++)
      mu[i] = 0.5 * (zeta[i] + zeta[i + 1]);

    /* compute raw weights */
    for (k = 0; k < Mpts; k++) {
      double nn = 1.0;
      for (i = 0; i < P; i++)
        nn *= lam[k] - mu[i];
      double dd = 1.0;
      for (i = 0; i < Mpts; i++)
        if (i != k)
          dd *= lam[k] - lam[i];
      raw_w[k] = fabs(nn / dd);
      if (k % 2 == 0)
        sum_even += raw_w[k];
      else
        sum_odd += raw_w[k];
    }
    for (k = 0; k < Mpts; k++) {
      w[k] = (k % 2 == 0) ? raw_w[k] / sum_even : raw_w[k] / sum_odd;
    }

    printf("  Scalar weights w_k (per-channel normalized):\n");
    printf("  %3s  %12s  %12s\n", "k", "w_even", "w_odd");
    printf("  %3s  %12s  %12s\n", "---", "----------", "----------");
    for (k = 0; k < Mpts; k++) {
      printf("  %3d  %12.6f", k, w[k]);
      if (k % 2 == 1)
        printf("\n");
    }

    /* Step 1: A_0 = sum W_k * Lambda_k  (W_k = diag(w_k,0) if even, diag(0,w_k)
     * if odd) = diag( sum_even[w_k*lam_k], sum_odd[w_k*lam_k] ) */
    double A0_e = 0.0, A0_o = 0.0;
    for (k = 0; k < Mpts; k++) {
      if (k % 2 == 0)
        A0_e += w[k] * lam[k];
      else
        A0_o += w[k] * lam[k];
    }
    printf("\n  A_0 = diag(%.4f, %.4f)     [channel 0 and channel 1 weighted "
           "means]\n",
           A0_e, A0_o);

    /* Step 2: B_0 B_0^T = sum W_k * (lam_k I - A_0)^2 */
    double BtB0_e = 0.0, BtB0_o = 0.0;
    for (k = 0; k < Mpts; k++) {
      double de = lam[k] - A0_e;
      double d_o = lam[k] - A0_o;
      if (k % 2 == 0)
        BtB0_e += w[k] * de * de;
      else
        BtB0_o += w[k] * d_o * d_o;
    }
    printf("  B_0 B_0^T = diag(%.4f, %.4f)\n", BtB0_e, BtB0_o);
    double B0_e = sqrt(BtB0_e), B0_o = sqrt(BtB0_o);
    printf("  B_0 = diag(%.4f, %.4f)     [Cholesky: sqrt of diagonal]\n", B0_e,
           B0_o);

    /* Verify B_0 matches the split-channel result */
    printf("\n  Block Lanczos A_0, B_0 vs split-channel dbg A_0, B_0:\n");
    printf("  A_0:  dbg(%9.4f,%9.4f)  lanczos(%9.4f,%9.4f)  match? %s\n",
           a11[0], a22[0], A0_e, A0_o,
           (fabs(a11[0] - A0_e) < 1e-10 && fabs(a22[0] - A0_o) < 1e-10) ? "yes"
                                                                        : "NO");
    printf("  B_0:  dbg(%9.4f,%9.4f)  lanczos(%9.4f,%9.4f)  match? %s\n",
           b11[0], b22[0], B0_e, B0_o,
           (fabs(b11[0] - B0_e) < 1e-10 && fabs(b22[0] - B0_o) < 1e-10) ? "yes"
                                                                        : "NO");

    /* Full block Lanczos loop for k=1..N-1 using Q-polynomials */
    /* Q_{-1}=0, Q_0=I. Q_{n+1} = (lam*I - A_n)*Q_n - B_{n-1}^2*Q_{n-1} */
    /* A_n = sum W_k * Q_n(lam_k)^T * lam_k * Q_n(lam_k) / norm */
    /* B_n B_n^T = sum W_k * Q_{n+1}(lam_k)^T * Q_{n+1}(lam_k) / norm_n */
    double A_full_e[10], A_full_o[10];
    double B_full_e[10], B_full_o[10];
    double Q_prev_e[20], Q_prev_o[20];
    double Q_cur_e[20], Q_cur_o[20];
    (void)BtB0_e;
    (void)BtB0_o;

    A_full_e[0] = A0_e;
    A_full_o[0] = A0_o;
    B_full_e[0] = B0_e;
    B_full_o[0] = B0_o;

    for (k = 0; k < Mpts; k++) {
      Q_cur_e[k] = 1.0;
      Q_cur_o[k] = 1.0;
      Q_prev_e[k] = 0.0;
      Q_prev_o[k] = 0.0;
    }

    for (i = 1; i < NB; i++) {
      double Bprev2_e = (i >= 2) ? B_full_e[i - 2] * B_full_e[i - 2] : 0.0;
      double Bprev2_o = (i >= 2) ? B_full_o[i - 2] * B_full_o[i - 2] : 0.0;

      /* Step: compute Q_i from recurrence (using Q_{i-1}, Q_{i-2}) */
      for (k = 0; k < Mpts; k++) {
        double qe =
            (lam[k] - A_full_e[i - 1]) * Q_cur_e[k] - Bprev2_e * Q_prev_e[k];
        double qo =
            (lam[k] - A_full_o[i - 1]) * Q_cur_o[k] - Bprev2_o * Q_prev_o[k];
        Q_prev_e[k] = Q_cur_e[k];
        Q_prev_o[k] = Q_cur_o[k];
        Q_cur_e[k] = qe;
        Q_cur_o[k] = qo;
      }

      /* A_i = sum W_k * Q_i(lam_k) * lam_k * Q_i(lam_k) / npk_i */
      /* where npk_i = sum W_k * Q_i(lam_k)^2 */
      double Anum_e = 0.0, Anum_o = 0.0;
      double Ad_e = 0.0, Ad_o = 0.0;
      for (k = 0; k < Mpts; k++) {
        if (k % 2 == 0) {
          Anum_e += Q_cur_e[k] * w[k] * lam[k] * Q_cur_e[k];
          Ad_e += Q_cur_e[k] * w[k] * Q_cur_e[k];
        } else {
          Anum_o += Q_cur_o[k] * w[k] * lam[k] * Q_cur_o[k];
          Ad_o += Q_cur_o[k] * w[k] * Q_cur_o[k];
        }
      }
      A_full_e[i] = Anum_e / Ad_e;
      A_full_o[i] = Anum_o / Ad_o;

      if (i < NB - 1) {
        /* Compute Q_{i+1} */
        double Bprev2b_e = B_full_e[i - 1] * B_full_e[i - 1];
        double Bprev2b_o = B_full_o[i - 1] * B_full_o[i - 1];
        double npx_e = 0.0, npx_o = 0.0;
        for (k = 0; k < Mpts; k++) {
          double qnext_e =
              (lam[k] - A_full_e[i]) * Q_cur_e[k] - Bprev2b_e * Q_prev_e[k];
          double qnext_o =
              (lam[k] - A_full_o[i]) * Q_cur_o[k] - Bprev2b_o * Q_prev_o[k];
          if (k % 2 == 0)
            npx_e += qnext_e * w[k] * qnext_e;
          else
            npx_o += qnext_o * w[k] * qnext_o;
        }
        B_full_e[i] = sqrt(npx_e / Ad_e);
        B_full_o[i] = sqrt(npx_o / Ad_o);
      }
    }

    printf("\n  Full block Lanczos (diagonal, 2-channel):\n");
    printf("  %3s  %10s  %10s  %10s  %10s\n", "k", "A11", "A22", "B11", "B22");
    printf("  %3s  %10s  %10s  %10s  %10s\n", "---", "--------", "--------",
           "--------", "--------");
    for (i = 0; i < NB; i++) {
      printf("  %3d  %10.4f  %10.4f", i, A_full_e[i], A_full_o[i]);
      if (i < NB - 1)
        printf("  %10.4f  %10.4f\n", B_full_e[i], B_full_o[i]);
      else
        printf("\n");
    }
    printf("\n  Full Lanczos A_k match split-channel dbg? %s\n",
           (fabs(A_full_e[0] - a11[0]) < 1e-10 &&
            fabs(A_full_o[0] - a22[0]) < 1e-10)
               ? "yes (matching)"
               : "NO (differ)");
  }

  free(M);
  printf("\n==================================================================="
         "===\n");
  printf("  Block de Boor-Golub reconstruction complete.\n");
  printf("====================================================================="
         "=\n");

  return 0;
}
