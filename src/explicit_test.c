/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Explicit formula test: spectral shift vs prime sum
 * @paper   yamaguchi-rh-2026.tex, Section 5.2
 * @theorem Theorem I (Guinand-Weil Explicit Formula)
 * @proof   N_J(E) - N_free(E) = (1/pi) sum sin(E log p)/sqrt(p)
 * @step    1 -- forward explicit formula verification
 */

/* explicit_test.c — Practical Application #5: Explicit Formula Test
 *
 * Tests the explicit formula:
 *   N_J(E) - N_free(E) ≈ (1/π) · Σ_{p ≤ Pmax} sin(E·log p)/√p
 *
 * Implementation:
 *   1. Reconstruct Jacobi matrix via de Boor-Golub with midpoints at N=25
 *   2. Compute eigenvalue counting function N_J(E) = #{k : γ_k ≤ E}
 *   3. Compute free counting function N_free(E) = #{k : g_k ≤ E} using Gram points
 *   4. Compute the spectral shift ξ(E) = N_J(E) - N_free(E)
 *   5. Compute ξ_pred(E) = (1/π) Σ_{p≤100} sin(E·log p)/√p
 *   6. Compute ξ(E) vs ξ_pred(E) at 20 energy points between γ_0 and γ_24
 *   7. Compute the correlation between ξ and ξ_pred
 *
 * Compile: gcc -O3 -fno-strict-aliasing -fno-peel-loops -fno-unswitch-loops
 *               -Wall -Wextra -Wconversion -Wshadow -Werror -Isrc -lm
 *               -o bin/explicit_test.exe src/explicit_test.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

#ifdef _WIN32
#include <windows.h>
#endif

/* ---- Zeta zeros (first 30) ---- */
static const double zeta_zeros[30] = {
    14.134725141734695, 21.022039638771556, 25.010857580145689,
    30.424876125859512, 32.935061587739192, 37.586178158825675,
    40.918719012147498, 43.327073280915002, 48.005150881167161,
    49.773832477672300, 52.970321477714464, 56.446247697063392,
    59.347044002602352, 60.831778524609810, 65.112544048081602,
    67.079810529494168, 69.546401711173985, 72.067157674481905,
    75.704690699083926, 77.144840068874799, 79.337375020249368,
    82.910380854086029, 84.735492980517051, 87.425274613125225,
    88.809111207634459, 92.491899270558491, 94.651344040519888,
    95.870634228245308, 98.831194218193687,101.317851005731384
};

/* ---- Riemann-Siegel theta function ---- */
static double theta_s(double t) {
    if (t <= 1.0) return -M_PI / 8.0;
    double x = t / (2.0 * M_PI), u = 1.0 / t, u2 = u * u, u4 = u2 * u2;
    return 0.5 * t * log(x) - 0.5 * t - M_PI / 8.0
           + u / 48.0
           + 7.0 * u * u2 / 5760.0
           + 31.0 * u4 / 80640.0
           + 127.0 * u4 * u2 / 430080.0
           + 2555.0 * u4 * u4 / 27525120.0
           + 1414477.0 * u4 * u4 * u2 / 18681062400.0;
}

static double theta_p(double t) {
    if (t <= 2.0 * M_PI) return 1.0;
    double t2 = t * t;
    return 0.5 * log(t / (2.0 * M_PI))
           - 1.0 / (24.0 * t2)
           + 7.0 / (960.0 * t2 * t2)
           + 31.0 / (8064.0 * t2 * t2 * t2);
}

/* ---- Gram point: solution of theta(g_n) = pi * n ---- */
static double gram_point(int n) {
    double g = (n == 0) ? 17.845599540410860
              : 2.0 * M_PI * (double)n / log((double)n + 1.0);
    for (int i = 0; i < 40; i++) {
        double f = theta_s(g) - M_PI * (double)n;
        double fp = theta_p(g);
        if (fabs(fp) < 1e-16) break;
        double dg = f / fp;
        g -= dg;
        if (fabs(dg) < 1e-15 * (1.0 + fabs(g))) break;
    }
    return g;
}

/* ---- Prime test ---- */
static int is_prime(int n) {
    if (n < 2) return 0;
    if (n == 2) return 1;
    if (n % 2 == 0) return 0;
    for (int d = 3; d * d <= n; d += 2)
        if (n % d == 0) return 0;
    return 1;
}

/* ---- de Boor-Golub: reconstruct Jacobi matrix from eigenvalues + midpoints ---- */
static int deboor(const double *lam, const double *mu, int N,
                  double *a, double *b) {
    for (int k = 0; k < N - 1; k++)
        if (!(lam[k] < mu[k] && mu[k] < lam[k + 1])) return -1;

    double w[30] = {0}, ws = 0;
    for (int k = 0; k < N; k++) {
        double nm = 1;
        for (int j = 0; j < N - 1; j++)
            nm *= lam[k] - mu[j];
        double dn = 1;
        for (int j = 0; j < N; j++)
            if (j != k) dn *= lam[k] - lam[j];
        w[k] = nm / dn;
        if (w[k] < 0) return -2;
        ws += w[k];
    }
    for (int k = 0; k < N; k++) w[k] /= ws;

    a[0] = 0;
    for (int i = 0; i < N; i++) a[0] += w[i] * lam[i];
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
                double b2 = (j > 0) ? b[j - 1] * b[j - 1] : 0;
                double pn = (lam[i] - a[j]) * pc - b2 * pp;
                pp = pc; pc = pn;
            }
            num += w[i] * lam[i] * pc * pc;
        }
        a[k] = num / npk;

        if (k < N - 1) {
            double npx = 0;
            for (int i = 0; i < N; i++) {
                double pp = 0, pc = 1;
                for (int j = 0; j <= k; j++) {
                    double b2 = (j > 0) ? b[j - 1] * b[j - 1] : 0;
                    double pn = (lam[i] - a[j]) * pc - b2 * pp;
                    pp = pc; pc = pn;
                }
                npx += w[i] * pc * pc;
            }
            b[k] = sqrt(npx / npk);
            npk = npx;
        }
    }
    return 0;
}

/* ---- Sturm counting: number of eigenvalues of Jacobi <= x ---- */
static int sturm_count(const double *d, const double *e, int N, double x) {
    int c = 0;
    double pp = 0.0, pc = 1.0;
    for (int k = 0; k < N; k++) {
        double ek = (k > 0) ? e[k - 1] : 0.0;
        double pn = (d[k] - x) * pc - ek * ek * pp;
        if (fabs(pn) > 1e100) pn = (pn > 0) ? 1e100 : -1e100;
        if (pc * pn < 0.0) c++;
        pp = pc; pc = pn;
    }
    return c;
}

/* ---- Jacobi rotation diagonalization ---- */
static void jacobi_diag(int n, const double *a, const double *b, double *ev) {
    double *M = calloc((size_t)(n * n), sizeof(double));
    for (int k = 0; k < n; k++) {
        M[k * n + k] = a[k];
        if (k < n - 1)
            M[k * n + k + 1] = M[(k + 1) * n + k] = b[k];
    }
    for (int sw = 0; sw < 80; sw++) {
        double moff = 0;
        for (int p = 0; p < n - 1; p++)
            for (int q = p + 1; q < n; q++)
                if (fabs(M[p * n + q]) > moff)
                    moff = fabs(M[p * n + q]);
        if (moff < 1e-15) break;
        for (int p = 0; p < n - 1; p++)
            for (int q = p + 1; q < n; q++) {
                double apq = M[p * n + q];
                double tol = 1e-16 * (fabs(M[p * n + p])
                           + fabs(M[q * n + q]) + 1);
                if (fabs(apq) < tol) continue;
                double app = M[p * n + p], aqq = M[q * n + q];
                double tau = (aqq - app) / (2 * apq);
                double t = (tau >= 0)
                    ? 1.0 / (tau + sqrt(1.0 + tau * tau))
                    : -1.0 / (-tau + sqrt(1.0 + tau * tau));
                double c = 1.0 / sqrt(1.0 + t * t), s = t * c;
                for (int i = 0; i < n; i++) {
                    double vip = M[i * n + p], viq = M[i * n + q];
                    M[i * n + p] = vip * c - viq * s;
                    M[i * n + q] = vip * s + viq * c;
                }
                for (int j = 0; j < n; j++) {
                    double vpj = M[p * n + j], vqj = M[q * n + j];
                    M[p * n + j] = vpj * c - vqj * s;
                    M[q * n + j] = vpj * s + vqj * c;
                }
            }
    }
    for (int k = 0; k < n; k++) ev[k] = M[k * n + k];
    for (int i = 1; i < n; i++) {
        double key = ev[i]; int j = i - 1;
        while (j >= 0 && ev[j] > key) { ev[j + 1] = ev[j]; j--; }
        ev[j + 1] = key;
    }
    free(M);
}

/* ---- Compute explicit formula prediction at energy E ---- */
static double xi_pred(double E, int Pmax) {
    double s = 0.0;
    for (int p = 2; p <= Pmax; p++) {
        if (!is_prime(p)) continue;
        s += sin(E * log((double)p)) / sqrt((double)p);
    }
    return s / M_PI;
}

/* ---- Pearson correlation coefficient ---- */
static double pearson_corr(const double *x, const double *y, int n) {
    double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
    for (int i = 0; i < n; i++) {
        sx += x[i]; sy += y[i];
        sxx += x[i] * x[i]; syy += y[i] * y[i];
        sxy += x[i] * y[i];
    }
    double num = (double)n * sxy - sx * sy;
    double den_x = (double)n * sxx - sx * sx;
    double den_y = (double)n * syy - sy * sy;
    if (den_x <= 0 || den_y <= 0) return 0.0;
    return num / sqrt(den_x * den_y);
}

/* ---- Spearman rank correlation ---- */
static void assign_ranks(const double *v, int n, double *ranks) {
    int *idx = malloc((size_t)n * sizeof(int));
    for (int i = 0; i < n; i++) idx[i] = i;
    double *sv = malloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) sv[i] = v[i];

    for (int i = 0; i < n - 1; i++)
        for (int j = i + 1; j < n; j++)
            if (sv[idx[i]] > sv[idx[j]]) {
                int t = idx[i]; idx[i] = idx[j]; idx[j] = t;
            }

    for (int i = 0; i < n; i++) {
        int j = i;
        while (j + 1 < n && sv[idx[j + 1]] == sv[idx[i]]) j++;
        double avg_rank = (double)(i + j) / 2.0 + 1.0;
        for (int k = i; k <= j; k++) ranks[idx[k]] = avg_rank;
        i = j;
    }
    free(idx); free(sv);
}

static double spearman_corr(const double *x, const double *y, int n) {
    double *rx = malloc((size_t)n * sizeof(double));
    double *ry = malloc((size_t)n * sizeof(double));
    assign_ranks(x, n, rx);
    assign_ranks(y, n, ry);
    double r = pearson_corr(rx, ry, n);
    free(rx); free(ry);
    return r;
}

/* ---- RMS error ---- */
static double rms_error(const double *a, const double *b, int n) {
    double ss = 0;
    for (int i = 0; i < n; i++) {
        double d = a[i] - b[i];
        ss += d * d;
    }
    return sqrt(ss / (double)n);
}

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

    int N = 25;
    int Pmax = 100;
    int n_points = 20;
    double lam[30], mu[29], a[30], b[29];

    for (int k = 0; k < N; k++) lam[k] = zeta_zeros[k];
    for (int k = 0; k < N - 1; k++)
        mu[k] = 0.5 * (zeta_zeros[k] + zeta_zeros[k + 1]);

    printf("======================================================================\n");
    printf("  PRACTICAL APPLICATION #5: EXPLICIT FORMULA TEST\n");
    printf("======================================================================\n\n");
    printf("  Hypothesis:\n");
    printf("    xi(E) = N_J(E) - N_free(E) ≈ (1/pi) * Σ_{p≤%d} sin(E·log p)/√p\n\n", Pmax);
    printf("  Where:\n");
    printf("    N_J(E)   = #{k : gamma_k ≤ E}  (zeta zero counting)\n");
    printf("    N_free(E) = #{k : g_k ≤ E}     (Gram point counting)\n");
    printf("    xi_pred(E) = (1/pi) * Σ_{p≤%d} sin(E·log p)/√p\n\n", Pmax);

    /* Step 1: de Boor-Golub reconstruction */
    printf("----------------------------------------------------------------------\n");
    printf("  STEP 1: Jacobi matrix reconstruction via de Boor-Golub (N=%d)\n", N);
    printf("----------------------------------------------------------------------\n\n");

    int ret = deboor(lam, mu, N, a, b);
    if (ret != 0) {
        fprintf(stderr, "ERROR: de Boor-Golub failed with code %d\n", ret);
        return 1;
    }

    printf("  Coefficients (a_k, b_k):\n");
    printf("  %3s  %12s  %12s\n", "k", "a_k", "b_k");
    printf("  %3s  %12s  %12s\n", "---", "---", "---");
    for (int k = 0; k < N; k++)
        printf("  %3d  %12.6f  %12.6f\n", k, a[k], (k < N - 1) ? b[k] : 0.0);

    /* Eigenvalue verification */
    double ev[30];
    jacobi_diag(N, a, b, ev);
    double max_err = 0;
    for (int k = 0; k < N; k++) {
        double err = fabs(ev[k] - lam[k]);
        if (err > max_err) max_err = err;
    }
    printf("\n  Max eigenvalue reconstruction error: %.1e\n", max_err);

    /* Step 2-3: Gram points */
    printf("\n----------------------------------------------------------------------\n");
    printf("  STEP 2-3: Gram points g_k for k = 0..%d\n", N - 1);
    printf("----------------------------------------------------------------------\n\n");
    printf("  %3s  %12s\n", "k", "g_k");
    printf("  %3s  %12s\n", "---", "---");
    double gram_vals[30];
    for (int k = 0; k < N; k++) {
        gram_vals[k] = gram_point(k);
        printf("  %3d  %12.6f\n", k, gram_vals[k]);
    }

    /* Step 4-5-6: Compute xi(E) and xi_pred(E) at 20 energy points */
    printf("\n----------------------------------------------------------------------\n");
    printf("  STEP 4-5-6: Spectral shift xi(E) vs xi_pred(E)\n");
    printf("              at %d points between gamma_0 and gamma_%d\n",
           n_points, N - 1);
    printf("----------------------------------------------------------------------\n\n");

    double E_lo = zeta_zeros[0] - 1.0;
    double E_hi = zeta_zeros[N - 1] + 1.0;
    double dE = (E_hi - E_lo) / (double)(n_points - 1);

    double xi_val[20], xi_pred_val[20], E_val[20];

    printf("  %4s  %10s  %10s  %10s  %10s  %10s\n",
           "idx", "E", "N_J(E)", "N_free(E)", "xi(E)", "xi_pred(E)");
    printf("  %4s  %10s  %10s  %10s  %10s  %10s\n",
           "---", "---", "---", "---", "---", "---");

    for (int i = 0; i < n_points; i++) {
        E_val[i] = E_lo + (double)i * dE;
        int n_j = sturm_count(a, b, N, E_val[i]);
        int n_free = 0;
        for (int k = 0; k < N; k++)
            if (gram_vals[k] <= E_val[i]) n_free++;
        xi_val[i] = (double)(n_j - n_free);
        xi_pred_val[i] = xi_pred(E_val[i], Pmax);
        printf("  %4d  %10.3f  %10d  %10d  %10.4f  %10.4f\n",
               i, E_val[i], n_j, n_free, xi_val[i], xi_pred_val[i]);
    }

    /* Step 7: Correlations */
    printf("\n======================================================================\n");
    printf("  STEP 7: CORRELATION ANALYSIS\n");
    printf("======================================================================\n\n");

    double r_pearson = pearson_corr(xi_val, xi_pred_val, n_points);
    double r_spearman = spearman_corr(xi_val, xi_pred_val, n_points);
    double rmse = rms_error(xi_val, xi_pred_val, n_points);

    double xi_mean = 0, xi_pred_mean = 0;
    for (int i = 0; i < n_points; i++) {
        xi_mean += xi_val[i];
        xi_pred_mean += xi_pred_val[i];
    }
    xi_mean /= (double)n_points;
    xi_pred_mean /= (double)n_points;

    double xi_std = 0, xi_pred_std = 0;
    for (int i = 0; i < n_points; i++) {
        double dx = xi_val[i] - xi_mean;
        double dy = xi_pred_val[i] - xi_pred_mean;
        xi_std += dx * dx;
        xi_pred_std += dy * dy;
    }
    xi_std = sqrt(xi_std / (double)n_points);
    xi_pred_std = sqrt(xi_pred_std / (double)n_points);

    printf("  Pearson correlation  r = %+.6f\n", r_pearson);
    printf("  Spearman correlation rho = %+.6f\n", r_spearman);
    printf("  RMS error                 = %.6f\n", rmse);
    printf("\n");
    printf("  xi(E) statistics:\n");
    printf("    mean = %+.4f,  std = %.4f\n", xi_mean, xi_std);
    printf("  xi_pred(E) statistics:\n");
    printf("    mean = %+.4f,  std = %.4f\n", xi_pred_mean, xi_pred_std);

    /* Detailed prime breakdown */
    printf("\n----------------------------------------------------------------------\n");
    printf("  Prime-by-prime contribution analysis\n");
    printf("  xi_pred(E) = (1/pi) * Σ_{p≤%d} sin(E·log p)/√p\n\n", Pmax);
    printf("  %5s  %10s  %14s\n", "p", "log p", "coeff (1/(pi√p))");
    printf("  %5s  %10s  %14s\n", "---", "---", "---");

    int n_primes = 0;
    double max_coeff = 0;
    int max_coeff_p = 0;
    for (int p = 2; p <= Pmax; p++) {
        if (!is_prime(p)) continue;
        n_primes++;
        double coeff = 1.0 / (M_PI * sqrt((double)p));
        if (coeff > max_coeff) { max_coeff = coeff; max_coeff_p = p; }
    }
    double sum_abs_coeff = 0;
    for (int p = 2; p <= Pmax; p++) {
        if (!is_prime(p)) continue;
        double coeff = 1.0 / (M_PI * sqrt((double)p));
        sum_abs_coeff += coeff;
        if (p <= 23 || p == max_coeff_p)
            printf("  %5d  %10.6f  %14.8f\n", p, log((double)p), coeff);
    }
    printf("  ...\n");
    printf("  %d primes ≤ %d, max coeff at p=%d: %.6f\n",
           n_primes, Pmax, max_coeff_p, max_coeff);
    printf("  Sum |coeff| = %.6f  (truncation error bound)\n", sum_abs_coeff);

    /* Dense sample correlation (100 points for robustness) */
    printf("\n----------------------------------------------------------------------\n");
    printf("  Dense sample check (100 points):\n");
    printf("----------------------------------------------------------------------\n");

    {
        int np = 100;
        double de2 = (E_hi - E_lo) / (double)(np - 1);
        double *xv = malloc((size_t)np * sizeof(double));
        double *yp = malloc((size_t)np * sizeof(double));
        for (int i = 0; i < np; i++) {
            double Ei = E_lo + (double)i * de2;
            int nj = sturm_count(a, b, N, Ei);
            int nf = 0;
            for (int k = 0; k < N; k++)
                if (gram_vals[k] <= Ei) nf++;
            xv[i] = (double)(nj - nf);
            yp[i] = xi_pred(Ei, Pmax);
        }
        double rp = pearson_corr(xv, yp, np);
        double rs = spearman_corr(xv, yp, np);
        double re = rms_error(xv, yp, np);
        printf("  Pearson (100 pts)  r = %+.6f\n", rp);
        printf("  Spearman (100 pts) rho = %+.6f\n", rs);
        printf("  RMS error  (100 pts)     = %.6f\n", re);
        free(xv); free(yp);
    }

    /* Final verdict */
    printf("\n======================================================================\n");
    printf("  RESULTS SUMMARY\n");
    printf("======================================================================\n\n");
    printf("  Explicit formula prediction vs reconstructed Jacobi spectral shift:\n");
    printf("\n");
    printf("    Pearson  r  = %+.4f  %s\n",
           r_pearson, fabs(r_pearson) > 0.5 ? "(strong)" :
           fabs(r_pearson) > 0.3 ? "(moderate)" : "(weak)");
    printf("    Spearman rho = %+.4f\n", r_spearman);
    printf("    RMS error    = %.4f\n", rmse);
    printf("\n");
    printf("  Interpretation:\n");
    printf("    The explicit formula provides an a priori prediction of the\n");
    printf("    spectral shift purely from prime numbers, without solving\n");
    printf("    any eigenvalue problem. The correlation between xi(E) and\n");
    printf("    xi_pred(E) measures how well prime arithmetic captures the\n");
    printf("    oscillatory structure of zeta zeros.\n");

    return 0;
}
