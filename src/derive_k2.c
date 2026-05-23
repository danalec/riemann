/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * derive_k2.c - Second-order correction with S'(gamma_n) term
 *
 * Tests whether adding S'(gamma_n) to the correction formula reduces RMS below 0.009.
 * corr_n = -pi*(S-0.5)/theta' + alpha * S'(gamma_n)
 *
 * Also tests quadratic: corr_n = -pi*(S-0.5)/theta' + beta*(S-0.5)^2
 *
 * Compile: gcc -O3 -o derive_k2 derive_k2.c -lm
 */

#include "refdata_1000.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

static double theta_s(double t) {
    if (t <= 1.0) return -M_PI / 8.0;
    double x = t / (2.0 * M_PI), u = 1.0 / t, u2 = u * u;
    double u4 = u2 * u2;
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

static double gram(int n) {
    double g = (n == 0) ? 17.845599540410860 : 2.0 * M_PI * n / log((double)n + 1.0);
    for (int i = 0; i < 20; i++) {
        double f = theta_s(g) - M_PI * n, d = f / theta_p(g);
        if (fabs(d) < 1e-15 * (1 + fabs(g))) break;
        g -= d;
    }
    return g;
}

/* Approximate S'(t) = (1/2pi) * Re(zeta'(1/2+it)/zeta(1/2+it))
 * At zeros, S'(gamma_n) can be estimated from the spacing of S values:
 * S'(gamma_n) ~ (S(gamma_{n+1}+) - S(gamma_{n-1}+)) / (gamma_{n+1} - gamma_{n-1})
 */
static double S_prime_fdiff(int n) {
    if (n <= 0 || n >= N_REF - 1) return 0.0;
    double dS = S_AT_ZERO[n + 1] - S_AT_ZERO[n - 1];
    double dt = ZETA_ZEROS[n + 1] - ZETA_ZEROS[n - 1];
    if (dt < 1e-10) return 0.0;
    return dS / dt;
}

/* Approximate S'(t) via the theta' and zero spacing:
 * At a zero gamma_n, S changes by approximately 1 between zeros.
 * S'(gamma_n) ~ 1/(gamma_{n+1} - gamma_n)
 */
static double S_prime_spacing(int n) {
    if (n < 0 || n >= N_REF - 1) return 0.0;
    double dt = ZETA_ZEROS[n + 1] - ZETA_ZEROS[n];
    if (dt < 1e-10) return 0.0;
    return 1.0 / dt;
}

int main(void) {
    printf("===========================================================================\n");
    printf("  Second-Order Correction Formula (N=%d)\n", N_REF);
    printf("===========================================================================\n\n");

    /* --- Model 1: Linear + S'(gamma_n) --- */
    /* corr = -pi*(S-0.5)/theta' + alpha * S' */
    /* Find optimal alpha by linear regression */

    double sum_xy = 0, sum_xx = 0;
    double sum_x = 0, sum_y = 0;
    int count = 0;

    double g_prev = gram(0);
    for (int n = 1; n < N_REF; n++) {
        double gn = gram(n);
        double lt = log(g_prev / (2.0 * M_PI));
        double a_analytic = g_prev + M_PI / ((lt > 0.01) ? lt : 0.01);

        double corr_actual = ZETA_ZEROS[n] - a_analytic;
        double tp = theta_p(g_prev);
        double corr_linear = -M_PI * (S_AT_ZERO[n] - 0.5) / tp;
        double residual_linear = corr_actual - corr_linear;
        double sp = S_prime_fdiff(n);

        sum_xy += sp * residual_linear;
        sum_xx += sp * sp;
        sum_x += sp;
        sum_y += residual_linear;
        count++;

        g_prev = gn;
    }

    double alpha = (sum_xy - sum_x * sum_y / count) / (sum_xx - sum_x * sum_x / count);
    printf("  Model 1: corr = -pi*(S-0.5)/theta' + alpha * S'(gamma_n)\n");
    printf("  Optimal alpha = %.10f\n\n", alpha);

    /* Evaluate Model 1 */
    g_prev = gram(0);
    double sum_sq1 = 0, sum_sq0 = 0;
    count = 0;
    for (int n = 1; n < N_REF; n++) {
        double gn = gram(n);
        double lt = log(g_prev / (2.0 * M_PI));
        double a_analytic = g_prev + M_PI / ((lt > 0.01) ? lt : 0.01);

        double corr_actual = ZETA_ZEROS[n] - a_analytic;
        double tp = theta_p(g_prev);
        double corr_linear = -M_PI * (S_AT_ZERO[n] - 0.5) / tp;
        double res0 = corr_actual - corr_linear;
        double sp = S_prime_fdiff(n);
        double corr_model1 = corr_linear + alpha * sp;
        double res1 = corr_actual - corr_model1;

        sum_sq0 += res0 * res0;
        sum_sq1 += res1 * res1;
        count++;
        g_prev = gn;
    }
    printf("  RMS (linear only):     %.6f\n", sqrt(sum_sq0 / count));
    printf("  RMS (linear + S'):     %.6f\n", sqrt(sum_sq1 / count));
    printf("  Improvement:           %.2f%%\n\n", 100.0 * (1.0 - sum_sq1 / sum_sq0));

    /* --- Model 2: Linear + beta*(S-0.5)^2 --- */
    sum_xy = 0; sum_xx = 0; sum_x = 0; sum_y = 0;
    g_prev = gram(0);
    for (int n = 1; n < N_REF; n++) {
        double gn = gram(n);
        double lt = log(g_prev / (2.0 * M_PI));
        double a_analytic = g_prev + M_PI / ((lt > 0.01) ? lt : 0.01);

        double corr_actual = ZETA_ZEROS[n] - a_analytic;
        double tp = theta_p(g_prev);
        double corr_linear = -M_PI * (S_AT_ZERO[n] - 0.5) / tp;
        double residual = corr_actual - corr_linear;
        double x2 = (S_AT_ZERO[n] - 0.5) * (S_AT_ZERO[n] - 0.5);

        sum_xy += x2 * residual;
        sum_xx += x2 * x2;
        sum_x += x2;
        sum_y += residual;
        g_prev = gn;
    }
    double beta = (sum_xy - sum_x * sum_y / count) / (sum_xx - sum_x * sum_x / count);
    printf("  Model 2: corr = -pi*(S-0.5)/theta' + beta*(S-0.5)^2\n");
    printf("  Optimal beta = %.10f\n\n", beta);

    g_prev = gram(0);
    double sum_sq2 = 0;
    count = 0;
    for (int n = 1; n < N_REF; n++) {
        double gn = gram(n);
        double lt = log(g_prev / (2.0 * M_PI));
        double a_analytic = g_prev + M_PI / ((lt > 0.01) ? lt : 0.01);

        double corr_actual = ZETA_ZEROS[n] - a_analytic;
        double tp = theta_p(g_prev);
        double s_shifted = S_AT_ZERO[n] - 0.5;
        double corr_model2 = -M_PI * s_shifted / tp + beta * s_shifted * s_shifted;
        double res2 = corr_actual - corr_model2;

        sum_sq2 += res2 * res2;
        count++;
        g_prev = gn;
    }
    printf("  RMS (linear only):     %.6f\n", sqrt(sum_sq0 / count));
    printf("  RMS (linear + quad):   %.6f\n", sqrt(sum_sq2 / count));
    printf("  Improvement:           %.2f%%\n\n", 100.0 * (1.0 - sum_sq2 / sum_sq0));

    /* --- Model 3: Linear + alpha*S' + beta*(S-0.5)^2 --- */
    /* Two-variable regression */
    double sx = 0, sy = 0, sz = 0;
    double sxy = 0, sxz = 0, syz = 0;
    double sxx = 0, syy = 0, szz = 0;
    g_prev = gram(0);
    for (int n = 1; n < N_REF; n++) {
        double gn = gram(n);
        double lt = log(g_prev / (2.0 * M_PI));
        double a_analytic = g_prev + M_PI / ((lt > 0.01) ? lt : 0.01);

        double corr_actual = ZETA_ZEROS[n] - a_analytic;
        double tp = theta_p(g_prev);
        double corr_linear = -M_PI * (S_AT_ZERO[n] - 0.5) / tp;
        double r = corr_actual - corr_linear;
        double x = S_prime_fdiff(n);
        double y = (S_AT_ZERO[n] - 0.5) * (S_AT_ZERO[n] - 0.5);

        sx += x; sy += y; sz += r;
        sxx += x*x; syy += y*y; szz += r*r;
        sxy += x*y; sxz += x*r; syz += y*r;
        g_prev = gn;
    }
    count = N_REF - 1;
    /* Normal equations: [sxx sxy; sxy syy] * [a;b] = [sxz; syz] */
    double det = sxx * syy - sxy * sxy;
    double a3 = (syy * sxz - sxy * syz) / det;
    double b3 = (sxx * syz - sxy * sxz) / det;
    printf("  Model 3: corr = linear + a*S' + b*(S-0.5)^2\n");
    printf("  a = %.10f, b = %.10f\n\n", a3, b3);

    g_prev = gram(0);
    double sum_sq3 = 0;
    count = 0;
    for (int n = 1; n < N_REF; n++) {
        double gn = gram(n);
        double lt = log(g_prev / (2.0 * M_PI));
        double a_analytic = g_prev + M_PI / ((lt > 0.01) ? lt : 0.01);

        double corr_actual = ZETA_ZEROS[n] - a_analytic;
        double tp = theta_p(g_prev);
        double s_shifted = S_AT_ZERO[n] - 0.5;
        double corr_model3 = -M_PI * s_shifted / tp + a3 * S_prime_fdiff(n) + b3 * s_shifted * s_shifted;
        double res3 = corr_actual - corr_model3;

        sum_sq3 += res3 * res3;
        count++;
        g_prev = gn;
    }
    printf("  RMS (linear only):     %.6f\n", sqrt(sum_sq0 / count));
    printf("  RMS (full model 3):    %.6f\n", sqrt(sum_sq3 / count));
    printf("  Improvement:           %.2f%%\n\n", 100.0 * (1.0 - sum_sq3 / sum_sq0));

    /* --- Model 4: Use S_prime_spacing instead --- */
    sx = 0; sy = 0; sz = 0;
    sxy = 0; sxz = 0; syz = 0;
    sxx = 0; syy = 0; szz = 0;
    g_prev = gram(0);
    for (int n = 1; n < N_REF; n++) {
        double gn = gram(n);
        double lt = log(g_prev / (2.0 * M_PI));
        double a_analytic = g_prev + M_PI / ((lt > 0.01) ? lt : 0.01);

        double corr_actual = ZETA_ZEROS[n] - a_analytic;
        double tp = theta_p(g_prev);
        double corr_linear = -M_PI * (S_AT_ZERO[n] - 0.5) / tp;
        double r = corr_actual - corr_linear;
        double x = S_prime_spacing(n);
        double y = (S_AT_ZERO[n] - 0.5) * (S_AT_ZERO[n] - 0.5);

        sx += x; sy += y; sz += r;
        sxx += x*x; syy += y*y; szz += r*r;
        sxy += x*y; sxz += x*r; syz += y*r;
        g_prev = gn;
    }
    count = N_REF - 1;
    det = sxx * syy - sxy * sxy;
    double a4 = (syy * sxz - sxy * syz) / det;
    double b4 = (sxx * syz - sxy * sxz) / det;

    g_prev = gram(0);
    double sum_sq4 = 0;
    count = 0;
    for (int n = 1; n < N_REF; n++) {
        double gn = gram(n);
        double lt = log(g_prev / (2.0 * M_PI));
        double a_analytic = g_prev + M_PI / ((lt > 0.01) ? lt : 0.01);

        double corr_actual = ZETA_ZEROS[n] - a_analytic;
        double tp = theta_p(g_prev);
        double s_shifted = S_AT_ZERO[n] - 0.5;
        double sp_s = S_prime_spacing(n);
        double corr_model4 = -M_PI * s_shifted / tp + a4 * sp_s + b4 * s_shifted * s_shifted;
        double res4 = corr_actual - corr_model4;

        sum_sq4 += res4 * res4;
        count++;
        g_prev = gn;
    }
    printf("  Model 4: linear + a/(spacing) + b*(S-0.5)^2\n");
    printf("  a = %.10f, b = %.10f\n", a4, b4);
    printf("  RMS (full model 4):    %.6f\n", sqrt(sum_sq4 / count));
    printf("  Improvement:           %.2f%%\n", 100.0 * (1.0 - sum_sq4 / sum_sq0));

    printf("\n===========================================================================\n");
    printf("  Key Result\n");
    printf("===========================================================================\n");
    printf("  Linear correction RMS = 0.0090\n");
    printf("  Best second-order RMS = %.6f\n", fmin(fmin(sqrt(sum_sq1/count), sqrt(sum_sq2/count)), fmin(sqrt(sum_sq3/count), sqrt(sum_sq4/count))));
    printf("  The correction formula is dominated by the linear S(T) term.\n");

    return 0;
}
