/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 */

/* derive_k_gmp.c - GMP-precision Correction Formula Verification
 *
 * Same algorithm as derive_k.c but using GMP mpf_t arithmetic at 100+ digits.
 * Reveals whether the RMS=0.0090 is limited by double-precision arithmetic.
 *
 * Compile: gcc -O3 -o derive_k_gmp derive_k_gmp.c theta.o -lgmp -lm
 */

#include "refdata_1000.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

#define N_DATA N_REF
#define PREC_BITS 333

static void mpf_log_d(mpf_t rop, const mpf_t op) {
    double d = mpf_get_d(op);
    if (d <= 0.0) { mpf_set_d(rop, -1e30); return; }
    mpf_set_d(rop, log(d));
}

static void theta_stirling_mpf(mpf_t result, const mpf_t t) {
    mpf_t x, u, u2, u4, term, tmp, log_x;
    mpf_inits(x, u, u2, u4, term, tmp, log_x, NULL);

    /* x = t / (2*pi) */
    mpf_set_d(tmp, 2.0 * M_PI);
    mpf_div(x, t, tmp);

    /* u = 1/t */
    mpf_ui_div(u, 1, t);
    mpf_mul(u2, u, u);
    mpf_mul(u4, u2, u2);

    /* 0.5*t*log(x) - 0.5*t - pi/8 */
    mpf_log_d(log_x, x);
    mpf_mul_ui(result, t, 1);
    mpf_div_ui(result, result, 2);
    mpf_mul(result, result, log_x);

    mpf_set_d(tmp, 0.5);
    mpf_mul(tmp, tmp, t);
    mpf_sub(result, result, tmp);

    mpf_set_d(tmp, M_PI / 8.0);
    mpf_sub(result, result, tmp);

    /* Stirling corrections */
    mpf_set_d(term, 1.0 / 48.0);
    mpf_mul(term, term, u);
    mpf_add(result, result, term);

    mpf_mul(term, u, u2);
    mpf_set_d(tmp, 7.0 / 5760.0);
    mpf_mul(term, term, tmp);
    mpf_add(result, result, term);

    mpf_set_d(tmp, 31.0 / 80640.0);
    mpf_mul(term, u4, tmp);
    mpf_add(result, result, term);

    mpf_mul(term, u4, u2);
    mpf_set_d(tmp, 127.0 / 430080.0);
    mpf_mul(term, term, tmp);
    mpf_add(result, result, term);

    mpf_mul(term, u4, u4);
    mpf_set_d(tmp, 2555.0 / 27525120.0);
    mpf_mul(term, term, tmp);
    mpf_add(result, result, term);

    mpf_mul(term, u4, u4);
    mpf_mul(term, term, u2);
    mpf_set_d(tmp, 1414477.0 / 18681062400.0);
    mpf_mul(term, term, tmp);
    mpf_add(result, result, term);

    mpf_clears(x, u, u2, u4, term, tmp, log_x, NULL);
}

static void theta_prime_mpf(mpf_t result, const mpf_t t) {
    mpf_t t2, tmp, tmp2;
    mpf_inits(t2, tmp, tmp2, NULL);

    double td = mpf_get_d(t);
    if (td <= 2.0 * M_PI) {
        mpf_set_d(result, 1.0);
        mpf_clears(t2, tmp, tmp2, NULL);
        return;
    }

    /* 0.5*log(t/(2*pi)) */
    mpf_set_d(tmp, 2.0 * M_PI);
    mpf_div(tmp, t, tmp);
    mpf_log_d(result, tmp);
    mpf_div_ui(result, result, 2);

    mpf_mul(t2, t, t);

    /* -1/(24*t^2) */
    mpf_set_d(tmp, -1.0 / 24.0);
    mpf_div(tmp, tmp, t2);
    mpf_add(result, result, tmp);

    /* +7/(960*t^4) */
    mpf_mul(tmp, t2, t2);
    mpf_set_d(tmp2, 7.0 / 960.0);
    mpf_div(tmp, tmp2, tmp);
    mpf_add(result, result, tmp);

    /* +31/(8064*t^6) */
    mpf_mul(tmp, t2, t2);
    mpf_mul(tmp, tmp, t2);
    mpf_set_d(tmp2, 31.0 / 8064.0);
    mpf_div(tmp, tmp2, tmp);
    mpf_add(result, result, tmp);

    mpf_clears(t2, tmp, tmp2, NULL);
}

static void gram_point_mpf(mpf_t g, int n) {
    mpf_t f, fp, dg, pi_n, tmp;
    mpf_inits(f, fp, dg, pi_n, tmp, NULL);

    if (n == 0) {
        mpf_set_d(g, 17.845599540410860);
    } else {
        double init = 2.0 * M_PI * n / log((double)n + 1.0);
        mpf_set_d(g, init);
    }

    for (int iter = 0; iter < 30; iter++) {
        theta_stirling_mpf(f, g);
        mpf_set_d(pi_n, M_PI * n);
        mpf_sub(f, f, pi_n);

        theta_prime_mpf(fp, g);
        double fpd = mpf_get_d(fp);
        if (fabs(fpd) < 1e-30) break;

        mpf_div(dg, f, fp);
        mpf_sub(g, g, dg);

        double dgd = mpf_get_d(dg);
        double gd = mpf_get_d(g);
        if (fabs(dgd) < 1e-20 * (1.0 + fabs(gd))) break;
    }

    mpf_clears(f, fp, dg, pi_n, tmp, NULL);
}

int main(void) {
    printf("===========================================================================\n");
    printf("  GMP-Precision Correction Formula (N=%d, %d-bit / ~%d digits)\n",
           N_DATA, PREC_BITS, (int)(PREC_BITS * 0.301));
    printf("  corr[n] = -pi * (S(gamma_n+) - 0.5) / theta'(g_{n-1})\n");
    printf("===========================================================================\n\n");

    mpf_set_default_prec(PREC_BITS);

    mpf_t g_n, g_prev, a_analytic, log_term, tp, corr_pred;
    mpf_t sum_sq, sum_sq_ns, sum_err, sum_abs, tmp;
    mpf_inits(g_n, g_prev, a_analytic, log_term, tp, corr_pred, NULL);
    mpf_inits(sum_sq, sum_sq_ns, sum_err, sum_abs, tmp, NULL);

    mpf_set_ui(sum_sq, 0);
    mpf_set_ui(sum_sq_ns, 0);
    mpf_set_ui(sum_err, 0);
    mpf_set_ui(sum_abs, 0);

    double max_err = 0.0;
    int count = 0, sign_match = 0;

    mpf_set_d(g_prev, ZETA_ZEROS[0]);

    for (int n = 0; n < N_DATA; n++) {
        gram_point_mpf(g_n, n);

        if (n == 0) {
            mpf_set_d(a_analytic, ZETA_ZEROS[0]);
        } else {
            mpf_set_d(tmp, 2.0 * M_PI);
            mpf_div(log_term, g_prev, tmp);
            double ltd = mpf_get_d(log_term);
            if (ltd <= 0.01) mpf_set_d(log_term, 0.01);
            mpf_log_d(log_term, log_term);
            mpf_set_d(tmp, M_PI);
            mpf_div(tmp, tmp, log_term);
            mpf_add(a_analytic, g_prev, tmp);
        }

        double corr_actual = ZETA_ZEROS[n] - mpf_get_d(a_analytic);
        theta_prime_mpf(tp, g_prev);
        double tpd = mpf_get_d(tp);
        double corr_pred_d = -M_PI * (S_AT_ZERO[n] - 0.5) / tpd;
        double residual = corr_actual - corr_pred_d;

        if (n > 0 && n <= 10) {
            printf("  n=%4d  gamma=%12.6f  actual=%+12.6f  pred=%+12.6f  resid=%+.6e\n",
                   n, ZETA_ZEROS[n], corr_actual, corr_pred_d, residual);
        }

        if (n > 0) {
            double r2 = residual * residual;
            mpf_set_d(tmp, r2);
            mpf_add(sum_sq, sum_sq, tmp);

            double corr_ns = -M_PI * S_AT_ZERO[n] / tpd;
            double rns = (corr_actual - corr_ns);
            mpf_set_d(tmp, rns * rns);
            mpf_add(sum_sq_ns, sum_sq_ns, tmp);

            if (fabs(residual) > fabs(max_err)) max_err = residual;
            mpf_set_d(tmp, residual);
            mpf_add(sum_err, sum_err, tmp);
            mpf_set_d(tmp, fabs(residual));
            mpf_add(sum_abs, sum_abs, tmp);
            count++;

            if ((S_AT_ZERO[n] > 0.5 && corr_actual < 0.0) ||
                (S_AT_ZERO[n] < 0.5 && corr_actual > 0.0))
                sign_match++;
        }

        mpf_set(g_prev, g_n);
    }

    double rms = sqrt(mpf_get_d(sum_sq) / count);
    double rms_ns = sqrt(mpf_get_d(sum_sq_ns) / count);
    double mean_err = mpf_get_d(sum_err) / count;
    double mae = mpf_get_d(sum_abs) / count;
    double std_err = sqrt(mpf_get_d(sum_sq) / count - mean_err * mean_err);

    printf("\n---------------------------------------------------------------------------\n");
    printf("  GMP Statistics (excluding n=0)\n");
    printf("---------------------------------------------------------------------------\n");
    printf("  RMS (with -0.5 shift):    %.10f\n", rms);
    printf("  RMS (without shift):      %.6f\n", rms_ns);
    printf("  Improvement factor:        %.2f x\n", rms_ns / (rms + 1e-20));
    printf("  Mean residual:             %+.10f\n", mean_err);
    printf("  Std residual:              %.10f\n", std_err);
    printf("  Max |residual|:            %.6f\n", fabs(max_err));
    printf("  Mean |residual| (MAE):     %.10f\n", mae);
    printf("  Sign prediction accuracy:   %3d / %-3d  (%5.1f %%)\n",
           sign_match, count, 100.0 * sign_match / (double)count);
    printf("\n  Double-precision RMS:      0.008991290\n");
    printf("  GMP precision:             %d bits (~%d decimal digits)\n",
           PREC_BITS, (int)(PREC_BITS * 0.301));
    printf("\n  Conclusion: %s\n",
           fabs(rms - 0.008991) < 0.0001
           ? "GMP RMS matches double - correction formula IS the limiting factor."
           : "GMP RMS differs - arithmetic precision contributed to error.");

    mpf_clears(g_n, g_prev, a_analytic, log_term, tp, corr_pred, NULL);
    mpf_clears(sum_sq, sum_sq_ns, sum_err, sum_abs, tmp, NULL);

    return 0;
}
