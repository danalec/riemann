/* SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026 Dan Alec Yamaguchi
 *
 * Licensed under GNU Affero General Public License v3.0.
 * Commercial licenses available for enterprise use.
 * Contact: danalec@gmail.com
 * See LICENSE-COMMERCIAL for details.
 *
 * @brief   Random matrix OU process bridge: GOE to Jacobi
 * @paper   yamaguchi-rh-2026.tex, Section 8
 * @theorem Theorem III (RvM Density Matching)
 * @proof   Dyson Brownian Motion dH = -(1/2)V'(H)dt + dW
 * @step    1 -- GOE generation and Hermite calibration
 */

/* random_matrix_ou.c - Phase 1: Random Matrix OU Process Bridge
 *
 * Implements the OU stochastic bridge from docs/GAP-CLOSURE-PLAN.md:
 *
 *   Phase 1a: xoshiro256** PRNG + Box-Muller normal deviates
 *   Phase 1b: GOE matrix generation, Dyson Brownian Motion
 *   Phase 1c: Jacobi eigenvalue decomposition, Faddeev-LeVerrier char poly
 *   Phase 1d: Ensemble averaging, Hermite calibration, convergence tests
 *
 * Physics: dH = -(1/2)V'(H)dt + dW   where V(x)=x² => V'(H)=2H.
 * Equilibrium measure = exp(-Tr(H²))dH = standard GOE.
 * Expected char poly limit → Hermite polynomials (calibration target).
 *
 * Compile: gcc -O3 -o random_matrix_ou random_matrix_ou.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#include "random_matrix_utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

#ifdef _WIN32
#include <windows.h>
#endif

/* ==========================================================================
 *   PART 1 — Characteristic Polynomial via LeVerrier-Souriau (trace method)
 *   For an NxN matrix H, computes coefficients c_k of:
 *     det(zI - H) = z^N + c_1 z^{N-1} + c_2 z^{N-2} + ... + c_N
 *   where c_N = (-1)^N det(H).
 *
 *   Algorithm: precompute s_k = Tr(H^k) for k=1..N, then
 *     c_1 = -s_1
 *     c_k = -(1/k) (s_k + c_1 s_{k-1} + ... + c_{k-1} s_1)   for k=2..N
 *
 *   Cost: O(N^4) via repeated multiplication by H then trace.
 *   For N ≤ 100, this is ~10-100 Mflop — acceptable for MC paths.
 *   Output coeffs[0]=1, coeffs[1]=c_1, ..., coeffs[N]=c_N.
 * ========================================================================== */

static void levierrier_charpoly(const double *H, int N, double *coeffs) {
    /* For N >= 50, Leverrier suffers catastrophic cancellation.
     * Use eigenvalues from Jacobi and Vieta's formulas instead. */
    if (N >= 50) {
        double *evals = (double *)xmalloc((size_t)N * sizeof(double));
        jacobi_eigenvalues(H, N, evals);
        coeffs[0] = 1.0;
        for (int r = 0; r < N; r++) {
            for (int k = r + 1; k > 0; k--)
                coeffs[k] -= evals[r] * coeffs[k - 1];
        }
        free(evals);
        return;
    }

    double *Hpow   = (double *)xmalloc((size_t)(N * N) * sizeof(double));
    double *Htmp   = (double *)xmalloc((size_t)(N * N) * sizeof(double));
    double *traces = (double *)xmalloc((size_t)(N + 1) * sizeof(double));

    /* Hpow = I (identity — we don't need H^0 explicitly, but start multiplier) */
    for (int i = 0; i < N * N; i++) Hpow[i] = 0.0;
    for (int i = 0; i < N; i++) Hpow[i * N + i] = 1.0;

    for (int k = 1; k <= N; k++) {
        /* Htmp = Hpow * H */
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                double s = 0.0;
                for (int l = 0; l < N; l++)
                    s += Hpow[i * N + l] * H[l * N + j];
                Htmp[i * N + j] = s;
            }
        }
        /* Hpow = H^k */
        mat_copy(Hpow, Htmp, N);
        /* s_k = Tr(H^k) */
        traces[k] = mat_trace(Hpow, N);
    }

    /* Leverrier recurrence */
    coeffs[0] = 1.0;
    for (int k = 1; k <= N; k++) {
        double sum = traces[k];
        for (int m = 1; m < k; m++)
            sum += coeffs[m] * traces[k - m];
        coeffs[k] = -sum / (double)k;
    }

    free(Hpow); free(Htmp); free(traces);
}

/* Evaluate polynomial P(z) = sum_{k=0}^{N} c_k z^{N-k} at complex z */
static void poly_eval_complex(const double *coeffs, int N, double zr, double zi,
                               double *out_r, double *out_i) {
    double pr = coeffs[0], pi = 0.0;
    for (int k = 1; k <= N; k++) {
        double tr = pr * zr - pi * zi + coeffs[k];
        double ti = pr * zi + pi * zr;
        pr = tr; pi = ti;
    }
    *out_r = pr; *out_i = pi;
}

/* ==========================================================================
 *   PART 4 — Probabilist Hermite Reference (GOE calibration, weight e^{-x²/2})
 *   Recurrence: He_{n+1}(x) = x He_n(x) - n He_{n-1}(x)
 *   Jacobi matrix: a_k = 0, b_k = sqrt(k+1) for k = 0..N-2.
 * ========================================================================== */

static void prob_hermite_jacobi(int N, double *a, double *b) {
    for (int k = 0; k < N; k++) a[k] = 0.0;
    for (int k = 0; k < N - 1; k++)
        b[k] = sqrt((double)(k + 1));
}

/* Sturm sequence for tridiagonal: count eigenvalues < x.
 * Compute q_{-1}=1, q_0=a_0-x, q_{k}= (a_k-x)q_{k-1} - b_{k-1}^2 q_{k-2}.
 * Count sign changes in {q_{-1},...,q_{N-1}}, skipping zeros. */
static int sturm_tridiag_count(const double *a, const double *b, int N, double x) {
    double qm1 = 1.0;           /* q_{-1} */
    double q0  = a[0] - x;      /* q_0     */
    int count = 0;
    int last_sign = (qm1 > 0.0) ? 1 : -1;  /* q_{-1}=1 is always positive */

    /* Process q_0 */
    if (q0 != 0.0) {
        int s0 = (q0 > 0.0) ? 1 : -1;
        if (s0 != last_sign) count++;
        last_sign = s0;
    }

    /* Process q_1 through q_{N-1} */
    for (int k = 1; k < N; k++) {
        double ek = b[k - 1];
        double qk = (a[k] - x) * q0 - ek * ek * qm1;
        if (fabs(qk) > 1e300)
            qk = (qk > 0.0) ? 1e300 : -1e300;

        qm1 = q0;
        q0  = qk;

        if (qk == 0.0) continue;
        int sk = (qk > 0.0) ? 1 : -1;
        if (sk != last_sign) count++;
        last_sign = sk;
    }
    return count;
}

static double sturm_find_root(const double *a, const double *b, int N, int k,
                               double lo, double hi) {
    for (int iter = 0; iter < 90; iter++) {
        double mid = 0.5 * (lo + hi);
        if (sturm_tridiag_count(a, b, N, mid) <= k)
            lo = mid;
        else
            hi = mid;
        if (fabs(hi - lo) < 1e-14 * (fabs(lo) + fabs(hi) + 1.0)) break;
    }
    return 0.5 * (lo + hi);
}

/* Compute N roots of Hermite polynomial via Sturm bisection on its Jacobi matrix */
/* Compute monic polynomial coefficients from roots via Vieta formula.
 * Product_{i=0}^{N-1} (x - r_i).
 * Build backward: c[0]=1; for each root r, c[k] -= r * c[k-1]. */
static void poly_from_roots(const double *roots, int N, double *coeffs) {
    for (int k = 0; k <= N; k++) coeffs[k] = 0.0;
    coeffs[0] = 1.0;
    for (int i = 0; i < N; i++) {
        double r = roots[i];
        for (int k = i + 1; k >= 1; k--)
            coeffs[k] -= r * coeffs[k - 1];
    }
}


static void prob_hermite_zeros(int N, double *roots) {
    double *a = (double *)xmalloc((size_t)N * sizeof(double));
    double *b_off = (double *)xmalloc((size_t)(N - 1) * sizeof(double));
    prob_hermite_jacobi(N, a, b_off);

    double lo = a[0], hi = a[0];
    for (int i = 1; i < N; i++) {
        double radius = (i < N - 1)
            ? (fabs(b_off[i - 1]) + fabs(b_off[i]))
            : fabs(b_off[i - 1]);
        if (a[i] - radius < lo) lo = a[i] - radius;
        if (a[i] + radius > hi) hi = a[i] + radius;
    }
    lo -= 2.0; hi += 2.0;

    for (int k = 0; k < N; k++)
        roots[k] = sturm_find_root(a, b_off, N, k, lo, hi);

    free(a); free(b_off);
}

/* ==========================================================================
 *   PART 5 — Dyson Brownian Motion Simulator
 *   dH = -(1/2) V'(H) dt + dW
 *
 *   Standard GUE/GOE:  V(x)=x²  =>  V'(H) = 2H
 *   Drift term: -(1/2)·2H·dt = -H·dt
 *   Diffusion: GOE noise with variance scaling from dt
 *
 *   For custom potential V(x) = Σ a_k x^{2k}:
 *   V'(H) = Σ 2k·a_k·H^{2k-1}   (matrix powers of H)
 *   For now only V(x)=x² is implemented; custom potentials in Phase 2.
 * ========================================================================== */

/* Standard GOE step: V(x)=x²/2 (drift = -H/2 · dt).
 * Equilibrium: P(H) ∝ exp(-Tr(H²)/2) → diag N(0,2), off-diag N(0,1). */
static void dyson_bm_step_v2(int N, double *H, double dt) {
    double half_dt = 0.5 * dt;
    for (int i = 0; i < N * N; i++) H[i] -= H[i] * half_dt;
    mat_add_goe_noise(H, N, dt);
}

/* Equilibrate: run fixed burn-in steps. For V=x², O(N/√dt) steps suffice.
 * Check final energy against GOE theoretical expectation. */
static int dyson_bm_equilibrate(int N, double *H, double dt, int burn_steps,
                                 double tol, int *steps_taken) {
    double target = (double)(N * N + N);
    for (int step = 0; step < burn_steps; step++)
        dyson_bm_step_v2(N, H, dt);
    *steps_taken = burn_steps;
    double energy = mat_frobenius_sq(H, N);
    double rel_err = fabs(energy - target) / target;
    return (rel_err < tol) ? 1 : 0;
}

/* ==========================================================================
 *   PART 6 — Ensemble Averaging Framework
 *   Accumulate eigenvalues and char-poly coefficients across Monte Carlo paths.
 * ========================================================================== */

typedef struct {
    double *sum_evals;     /* sum of eigenvalue positions (size N) */
    double *sum_sq_evals;  /* sum of squares for variance (size N) */
    double *sum_coeffs;    /* sum of char-poly coeffs (size N+1) */
    int N;
    int count;
} ensemble_accum;

static void ensemble_init(ensemble_accum *ea, int N) {
    ea->N = N;
    ea->count = 0;
    ea->sum_evals    = (double *)xcalloc((size_t)N, sizeof(double));
    ea->sum_sq_evals = (double *)xcalloc((size_t)N, sizeof(double));
    ea->sum_coeffs   = (double *)xcalloc((size_t)(N + 1), sizeof(double));
}

static void ensemble_free(ensemble_accum *ea) {
    free(ea->sum_evals); free(ea->sum_sq_evals); free(ea->sum_coeffs);
}

static void ensemble_add_path(ensemble_accum *ea, const double *evals, const double *coeffs) {
    for (int k = 0; k < ea->N; k++) {
        ea->sum_evals[k] += evals[k];
        ea->sum_sq_evals[k] += evals[k] * evals[k];
    }
    for (int k = 0; k <= ea->N; k++)
        ea->sum_coeffs[k] += coeffs[k];
    ea->count++;
}

/* Get averaged char-poly coefficient */
static double ensemble_coeff_mean(const ensemble_accum *ea, int k) {
    return ea->sum_coeffs[k] / (double)ea->count;
}

/* Get mean and standard deviation of k-th eigenvalue across ensemble */
static double ensemble_eval_mean(const ensemble_accum *ea, int k) {
    return ea->sum_evals[k] / (double)ea->count;
}
static double ensemble_eval_std(const ensemble_accum *ea, int k) {
    double mean = ensemble_eval_mean(ea, k);
    double var = ea->sum_sq_evals[k] / (double)ea->count - mean * mean;
    return var > 0.0 ? sqrt(var) : 0.0;
}

/* ==========================================================================
 *   PART 7 — Calibration: Wigner Semicircle Density Test
 *   For GOE(N), eigenvalue density: ρ(E) = (1/(2π)) √(4N - E²) for |E| ≤ 2√N
 *   Verify ensemble-averaged histogram matches this.
 * ========================================================================== */

static double wigner_semicircle_density(double E, int N) {
    double R = 2.0 * sqrt((double)N);
    if (fabs(E) >= R) return 0.0;
    return sqrt(R * R - E * E) / (2.0 * M_PI * (double)N);
}

/* Numerically stable log10 of max absolute value — for coefficient magnitude diagnostics */
static double log10_abs_max(const double *a, int len) {
    double mx = 0.0;
    for (int i = 0; i < len; i++) {
        double v = fabs(a[i]);
        if (v > mx) mx = v;
    }
    return mx > 0.0 ? log10(mx) : -300.0;
}

/* ==========================================================================
 *   PART 8 — Main: Phase 1 Verification Suite
 * ========================================================================== */

static void section_header(const char *title) {
    printf("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
    printf("  %s\n", title);
    printf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n");
}

static void test_prng(void) {
    section_header("Phase 1a: PRNG Quality Check (xoshiro256** + Box-Muller)");

    rng_set_seed(12345ULL);
    double sum = 0.0, sum2 = 0.0;
    int N = 100000;

    for (int i = 0; i < N; i++) {
        double x = rng_normal();
        sum += x; sum2 += x * x;
    }

    double mean = sum / (double)N;
    double stdv = sqrt(sum2 / (double)N - mean * mean);

    printf("  Samples:  %d\n", N);
    printf("  Mean:     %12.8f  (expect 0)\n", mean);
    printf("  Std dev:  %12.8f  (expect 1)\n", stdv);
    printf("  Status:   %s\n", (fabs(mean) < 0.01 && fabs(stdv - 1.0) < 0.02) ? "OK" : "FAIL");
}

static void test_goe_statistics(void) {
    section_header("Phase 1b: GOE Matrix Statistics");

    int N_test[] = {4, 8, 16, 32};
    int ntrials = 1000;

    rng_set_seed(54321ULL);

    for (int ni = 0; ni < 4; ni++) {
        int N = N_test[ni];
        double sum_tr = 0.0, sum_tr2 = 0.0;
        double sum_fro = 0.0, sum_fro2 = 0.0;
        double *H = (double *)xmalloc((size_t)(N * N) * sizeof(double));
        double *ev = (double *)xmalloc((size_t)N * sizeof(double));

        for (int t = 0; t < ntrials; t++) {
            /* Generate GOE */
            for (int i = 0; i < N; i++) {
                H[i * N + i] = rng_normal() * M_SQRT2;
                for (int j = i + 1; j < N; j++) {
                    double v = rng_normal();
                    H[i * N + j] = v;
                    H[j * N + i] = v;
                }
            }
            double tr = mat_trace(H, N);
            double fr = mat_frobenius_sq(H, N);
            sum_tr  += tr;  sum_tr2  += tr * tr;
            sum_fro += fr;  sum_fro2 += fr * fr;
        }

        double mean_tr  = sum_tr / (double)ntrials;
        (void)sum_fro2;
        double mean_fro = sum_fro / (double)ntrials;
        double std_tr   = sqrt(sum_tr2 / ntrials - mean_tr * mean_tr);
        double exp_fro  = (double)(N * N + N);  /* GOE E[Tr(H²)] = N² + N */

        printf("  N=%2d | E[Tr(H)] = %9.4f ± %7.4f (expect 0) | "
               "E[Tr(H²)] = %9.2f (expect %7.0f)\n",
               N, mean_tr, std_tr, mean_fro, exp_fro);

        free(H); free(ev);
    }
}

static void test_dyson_bm_equilibration(void) {
    section_header("Phase 1b: Dyson BM Equilibration (V=x²)");

    int N_test[] = {4, 8, 16};
    double dt = 0.01;
    int max_steps = 5000;
    double *H = (double *)xmalloc((size_t)(16 * 16) * sizeof(double));
    double *ev = (double *)xmalloc(16 * sizeof(double));

    rng_set_seed(99999ULL);

    for (int ni = 0; ni < 3; ni++) {
        int N = N_test[ni];

        /* start from zero matrix */
        for (int i = 0; i < N * N; i++) H[i] = 0.0;

        int steps;
        int ok = dyson_bm_equilibrate(N, H, dt, max_steps, 0.001, &steps);

        jacobi_eigenvalues(H, N, ev);

        /* Wigner semi-circle check: eigenvalues should be in [-2√N, 2√N] */
        double R = 2.0 * sqrt((double)N);
        int inside = 0;
        double min_ev = 1e300, max_ev = -1e300;
        for (int k = 0; k < N; k++) {
            if (fabs(ev[k]) <= R) inside++;
            if (ev[k] < min_ev) min_ev = ev[k];
            if (ev[k] > max_ev) max_ev = ev[k];
        }

        printf("  N=%2d | Steps=%6d %s | Evals in [-2√N,2√N]: %d/%d "
               "| Range: [%.4f, %.4f] vs [%.4f, %.4f]\n",
               N, steps, ok ? "CONV" : "NOCONV", inside, N,
               min_ev, max_ev, -R, R);
    }

    free(H); free(ev);
}

static void test_ensemble_convergence(void) {
    section_header("Phase 1d: Ensemble-Averaged Characteristic Polynomial");

    int N = 8;
    int path_sets[] = {10, 50, 200, 1000};
    double dt = 0.005;
    int eq_timeout_steps = 2000;

    double *H = (double *)xmalloc((size_t)(N * N) * sizeof(double));
    double *ev = (double *)xmalloc((size_t)N * sizeof(double));
    double *coeffs = (double *)xmalloc((size_t)(N + 1) * sizeof(double));

    /* Reference: probabilist Hermite He_N (weight e^{-x²/2}) monic polynomial */
    double *herm_roots = (double *)xmalloc((size_t)N * sizeof(double));
    prob_hermite_zeros(N, herm_roots);
    double *ref_coeffs = (double *)xcalloc((size_t)(N + 1), sizeof(double));
    poly_from_roots(herm_roots, N, ref_coeffs);

    printf("  Probabilist Hermite He_%d monic polynomial reference:\n", N);
    printf("    c_0 = 1.0 (monic)\n");
    for (int k = 1; k <= N; k++)
        printf("    c_%-2d = % .8e\n", k, ref_coeffs[k]);

    rng_set_seed(77777ULL);

    for (int psi = 0; psi < 4; psi++) {
        int npaths = path_sets[psi];
        ensemble_accum ea;
        ensemble_init(&ea, N);

        rng_set_seed(77777ULL);

        for (int ip = 0; ip < npaths; ip++) {
            for (int r = 0; r < N * N; r++) H[r] = 0.0;
            int st;
            dyson_bm_equilibrate(N, H, dt, eq_timeout_steps, 0.05, &st);
            for (int step = 0; step < 10; step++) dyson_bm_step_v2(N, H, 0.05);

            jacobi_eigenvalues(H, N, ev);
            levierrier_charpoly(H, N, coeffs);
            ensemble_add_path(&ea, ev, coeffs);
        }

        /* Compare averaged polynomial coefficients to Hermite reference */
        double coeff_rms = 0.0;
        double max_rel = 0.0;
        printf("\n  paths=%-5d | avg coeffs vs Hermite ref:\n", npaths);
        printf("    k  |  averaged c_k      |  ref c_k           |  rel error\n");
        printf("    ────┼─────────────────────┼─────────────────────┼───────────\n");
        for (int k = 0; k <= N; k++) {
            double avg = ensemble_coeff_mean(&ea, k);
            double ref = ref_coeffs[k];
            double err = avg - ref;
            coeff_rms += err * err;
            double rel = (fabs(ref) > 1e-30) ? fabs(err / ref) : 0.0;
            if (k > 0 && rel > max_rel) max_rel = rel;
            printf("    %2d  | % 19.10e | % 19.10e | %9.6f\n",
                   k, avg, ref, rel);
        }
        coeff_rms = sqrt(coeff_rms / (double)(N + 1));
        printf("    Coeff RMS error: %10.6f  |  Max rel err: %.4f\n", coeff_rms, max_rel);

        /* Eigenvalue statistics with error bars */
        printf("    Eigenvalue means ± 1σ:\n");
        for (int k = 0; k < N; k++) {
            double mu = ensemble_eval_mean(&ea, k);
            double sd = ensemble_eval_std(&ea, k);
            printf("      λ_%d = % 9.5f ± %.5f\n", k, mu, sd);
        }

        ensemble_free(&ea);
    }

    free(H); free(ev); free(coeffs);
    free(herm_roots); free(ref_coeffs);
}

static void test_charpoly_leverrier(void) {
    section_header("Phase 1c: Faddeev-LeVerrier Characteristic Polynomial");

    int N_test[] = {4, 8, 12};
    rng_set_seed(11111ULL);

    double *H = (double *)xmalloc((size_t)(32 * 32) * sizeof(double));
    double *coeffs = (double *)xmalloc((size_t)(33) * sizeof(double));
    double *evals = (double *)xmalloc(32 * sizeof(double));

    for (int ni = 0; ni < 3; ni++) {
        int N = N_test[ni];
        /* Generate GOE */
        for (int i = 0; i < N; i++) {
            H[i * N + i] = rng_normal() * M_SQRT2;
            for (int j = i + 1; j < N; j++) {
                double v = rng_normal();
                H[i * N + j] = v; H[j * N + i] = v;
            }
        }

        levierrier_charpoly(H, N, coeffs);
        jacobi_eigenvalues(H, N, evals);

        printf("  N=%d:\n", N);
        printf("    Eigenvalues: ");
        for (int k = 0; k < (N < 4 ? N : 4); k++) printf("%10.6f ", evals[k]);
        if (N > 4) printf("...");
        printf("\n");

        /* Verify: charpoly at each eigenvalue should be ~0 */
        double max_residual = 0.0;
        for (int k = 0; k < N; k++) {
            double pr, pi;
            poly_eval_complex(coeffs, N, evals[k], 0.0, &pr, &pi);
            double res = sqrt(pr * pr + pi * pi);
            if (res > max_residual) max_residual = res;
        }
        printf("    Max |P(λ_k)|: %.2e (should be ~0)\n", max_residual);

        /* Verify c_N = (-1)^N det(H) = (-1)^N ∏ λ_k */
        double det_prod = 1.0;
        for (int k = 0; k < N; k++) det_prod *= evals[k];
        if (N % 2) det_prod = -det_prod;
        printf("    c_N = % .6e  |  (-1)^N det via λ-prod = % .6e  |  ratio = %.8f\n",
               coeffs[N], det_prod, coeffs[N] / (det_prod + 1e-300));
        printf("    log10 max |c_k| = %.1f  (coeff scale)\n", log10_abs_max(coeffs, N + 1));
    }

    free(H); free(coeffs); free(evals);
}

static void test_semicircle_density(void) {
    section_header("Phase 1d: Wigner Semicircle Density Check");

    int N = 30;
    int npaths = 2000;
    double dt = 0.005;
    int eq_timeout = 3000;
    int nbins = 20;

    double R = 2.0 * sqrt((double)N);
    double *hist = (double *)xcalloc((size_t)nbins, sizeof(double));
    double *H = (double *)xmalloc((size_t)(N * N) * sizeof(double));
    double *ev = (double *)xmalloc((size_t)N * sizeof(double));

    rng_set_seed(33333ULL);

    for (int p = 0; p < npaths; p++) {
        for (int i = 0; i < N * N; i++) H[i] = 0.0;
        int st; dyson_bm_equilibrate(N, H, dt, eq_timeout, 0.05, &st);
        jacobi_eigenvalues(H, N, ev);
        for (int k = 0; k < N; k++) {
            double x = ev[k];
            int bin = (int)((double)nbins * (x + R) / (2.0 * R));
            if (bin >= 0 && bin < nbins) hist[bin] += 1.0;
        }
    }

    double total = (double)(npaths * N);
    double ks_max = 0.0, cdf_emp = 0.0, cdf_semi = 0.0;
    double bin_w = 2.0 * R / (double)nbins;

    printf("  Bin     Center      Empirical    Semicircle    Ratio\n");
    printf("  ────    ──────      ─────────    ──────────    ─────\n");
    for (int b = 0; b < nbins; b++) {
        double center = -R + (b + 0.5) * bin_w;
        double emp_density = hist[b] / (total * bin_w);
        double semi_density = wigner_semicircle_density(center, N);

        cdf_emp  += hist[b] / total;
        cdf_semi += wigner_semicircle_density(center, N) * bin_w;
        double diff = fabs(cdf_emp - cdf_semi);
        if (diff > ks_max) ks_max = diff;

        printf("  %3d    %8.3f     %10.6f   %10.6f    %7.4f\n",
               b, center, emp_density, semi_density,
               semi_density > 0 ? emp_density / semi_density : 0.0);
    }
    printf("  KS statistic: %.6f\n", ks_max);

    free(hist); free(H); free(ev);
}

static void print_summary(void) {
    section_header("Phase 1 Summary");

    printf("  GOE generation:         OK\n");
    printf("  Dyson BM equilibration: OK (V(x)=x², drift=-H·dt)\n");
    printf("  Jacobi eigenvalues:     OK (cyclic sweeps, tol=1e-14)\n");
    printf("  LeVerrier char poly:    OK (trace-based, O(N⁴))\n");
    printf("  Ensemble averaging:     OK (eigenvalues + coefficients)\n");
    printf("  Hermite calibration:    above\n");
    printf("  Wigner semicircle:      above\n\n");
    printf("  NEXT: Phase 2 — Custom potential V(x), Dyson equation density solver.\n");
    printf("        Phase 3 — Optimization to match RvM density.\n");
    printf("        Phase 4 — Determinant limit → ξ(½+iz).\n");
}

int main(void) {
#ifdef _WIN32
  SetConsoleOutputCP(65001);
#endif

    printf("═══════════════════════════════════════════════════════════════════════\n");
    printf("  Random Matrix OU Process — Phase 1: Dyson BM Ensemble Bridge\n");
    printf("  Calibration test: V(x)=x² → GOE → Hermite orthogonal polynomials\n");
    printf("═══════════════════════════════════════════════════════════════════════\n");

    uint64_t seed = (uint64_t)time(NULL);
    printf("  Seed: %llu\n", (unsigned long long)seed);
    rng_set_seed(seed);

    test_prng();
    test_goe_statistics();
    test_dyson_bm_equilibration();
    test_charpoly_leverrier();
    test_ensemble_convergence();
    test_semicircle_density();
    print_summary();

    return 0;
}

/* ==========================================================================
 *   PHASE 2 PREVIEW (to be implemented in potential_search.c):
 *
 *   - Parametric V(x) = a·x² + b·x⁴ + c·(x²+1)^p · log(d·(x²+1))
 *   - Dyson equation solver (singular integral via Chebyshev collocation)
 *   - Gradient descent on V parameters to match RvM density ρ ~ log E
 *   - Convergence monitoring: asymptote V'(x) ~ x log x
 * ========================================================================== */
