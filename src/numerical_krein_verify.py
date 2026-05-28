# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dan Alec Yamaguchi
#
# Licensed under GNU Affero General Public License v3.0.
# Commercial licenses available for enterprise use.
# Contact: danalec@gmail.com
# See LICENSE-COMMERCIAL for details.
#
"""
Numerical verification of spectral shift identity and eigenvalue-zero correspondence.

Key comparisons:
1. Eigenvalues of J_N vs known Riemann zeta zeros _k
2. Spectral shift function: is Tr h(J_N) - Tr h(J_N^free) approaching the GW prime sum?
3. Per-prime decomposition of eigenvalue shifts
4. Krein determinant arg det(I + G P) vs spectral shift function
5. Off-diagonal perturbation contribution (unexplored direction)
"""

import numpy as np
from scipy import special
from scipy.linalg import eigvalsh_tridiagonal
from sympy import primerange
import warnings
warnings.filterwarnings('ignore')

def theta(t):
    s = 0.25 + 0.5j * t
    return np.imag(special.loggamma(s)) - 0.5 * t * np.log(np.pi)

def theta_prime(t):
    return 0.5 * np.log(t / (2 * np.pi))

def gram_point(n):
    if n == 0:
        t = 17.85
    elif n == 1:
        t = 23.0
    elif n == 2:
        t = 28.0
    else:
        t = 2 * np.pi * n / np.log(n)
    for _ in range(200):
        th = theta(t)
        thp = theta_prime(t)
        if abs(thp) < 1e-15:
            break
        dt = (n * np.pi - th) / thp
        t += dt
        if t <= 0:
            t = 1.0
        if abs(dt) < 1e-13:
            break
    return t

# First 80 known Riemann zeta zeros (imaginary parts)
ZETA_ZEROS = [
    14.134725141734693, 21.022039638771554, 25.010857580145688,
    30.424876125859513, 32.935061587739190, 37.586178158825671,
    40.918719012147495, 43.327073280914999, 48.005150881167159,
    49.773832477672302, 52.970321477714460, 56.446247697063395,
    59.347044002602903, 60.831778524609809, 65.112544048081607,
    67.079810529501173, 69.546401711173979, 72.103742237830137,
    75.704690699083933, 77.144840068874805, 79.337375020923417,
    82.910380654236103, 84.735492980517051, 87.425274613125220,
    88.809111207634466, 92.491899270558484, 94.651344025671971,
    95.870634228245310, 98.831194218193692, 101.317851005731391,
    103.725538040478337, 105.446623052326849, 107.168611184276407,
    111.029535543169681, 111.874659177886826, 114.320220914455177,
    116.226696320656212, 118.790766072935096, 121.370125002420689,
    122.946829293552122, 124.256818924159199, 127.516683879596500,
    129.578704199956060, 131.087688530932655, 133.497737202997583,
    134.756509753000232, 138.116042054676727, 139.736208952084367,
    141.123707404050607, 143.111845807612696, 146.000982486765354,
    147.422765342557332, 150.053520420545660, 150.925257612771248,
    153.024693081921520, 155.810957307765530, 158.076777823405904,
    159.092586242846180, 161.188961138241624, 163.030709687267800,
    165.537069335201030, 167.184439715425005, 169.094515448073979,
    172.136843090947076, 174.754191824854660, 176.441434297749416,
    178.377050085186995, 179.916484021710359, 182.207028392384838,
    184.877464796349377, 185.598783020854500, 187.228933012046095,
    189.416779930056275, 192.026656360710506, 193.079726801190097,
    195.265396939468963, 196.876481830399611, 198.015309672474026,
    201.264751541703881, 202.493594516780866, 204.189671109100816,
    207.124903734991267, 209.117306730210095,
]

def build_free_jacobi(N):
    g = np.array([gram_point(n) for n in range(N + 1)])
    a = np.zeros(N)
    b = np.zeros(N - 1)
    for n in range(N):
        a[n] = g[n - 1] + np.pi / np.log(g[n - 1] / (2 * np.pi)) if n > 0 else g[0] + np.pi / np.log(g[0] / (2 * np.pi))
    for n in range(N - 1):
        spacing = g[n + 1] - g[n]
        dens = theta_prime(g[n + 1])
        b[n] = np.sqrt(spacing) * dens
    return a, b, g

def build_single_prime_perturbation(N, g, p):
    alpha_p = -np.log(p) / (2 * np.pi * np.sqrt(p))
    logp = np.log(p)
    dP = np.zeros(N)
    oP = np.zeros(N - 1)
    for n in range(N):
        dP[n] = alpha_p * np.sin(logp * g[n])
    for n in range(N - 1):
        g_avg = (g[n] + g[n + 1]) / 2
        oP[n] = alpha_p * np.cos(logp * g_avg)
    return dP, oP, alpha_p

def build_full_perturbation(N, g, Lambda):
    dP = np.zeros(N)
    oP = np.zeros(N - 1)
    for p in primerange(2, int(Lambda) + 1):
        dP_p, oP_p, _ = build_single_prime_perturbation(N, g, p)
        dP += dP_p
        oP += oP_p
    return dP, oP

def main():
    print("=" * 80)
    print("NUMERICAL SPECTRAL SHIFT VERIFICATION")
    print("=" * 80)
    
    for N in [30, 50, 80]:
        print(f"\n{'='*80}")
        print(f"N = {N}")
        print(f"{'='*80}")
        
        a_free, b_free, g = build_free_jacobi(N)
        Lambda = g[N - 1]
        
        print(f"Gram points: g_0={g[0]:.4f}, g_{{N-1}}={g[N-1]:.4f}")
        print(f"Lambda = {Lambda:.2f}, #primes <= Lambda = {len(list(primerange(2, int(Lambda)+1)))}")
        
        eigs_free = eigvalsh_tridiagonal(a_free, b_free)
        
        dP, oP = build_full_perturbation(N, g, Lambda)
        
        gamma_1 = ZETA_ZEROS[0]
        Delta = gamma_1 - (g[0] + np.pi / np.log(g[0] / (2 * np.pi)))
        
        a_full = a_free.copy()
        a_full[0] += Delta
        a_full += dP
        b_full = b_free + oP
        
        eigs_full = eigvalsh_tridiagonal(a_full, b_full)
        
        # ========================================
        # KEY COMPARISON 1: Eigenvalues vs Zeta Zeros
        # ========================================
        n_compare = min(N, len(ZETA_ZEROS))
        print(f"\n--- Eigenvalue-Zero Comparison (first {n_compare}) ---")
        print(f"{'k':>4} {'gamma_k':>14} {'lambda_k':>14} {'lambda_free':>14} {'err(full)':>12} {'err(free)':>12}")
        
        rms_full = 0.0
        rms_free = 0.0
        count = 0
        for k in range(n_compare):
            gamma_k = ZETA_ZEROS[k]
            lam_k = eigs_full[k]
            lam_free_k = eigs_free[k]
            err_full = lam_k - gamma_k
            err_free = lam_free_k - gamma_k
            if abs(gamma_k) > 1e-10:
                rms_full += err_full**2
                rms_free += err_free**2
                count += 1
            if k < 15 or k >= n_compare - 5:
                print(f"{k:4d} {gamma_k:14.6f} {lam_k:14.6f} {lam_free_k:14.6f} {err_full:12.6f} {err_free:12.6f}")
            elif k == 15:
                print("  ...")
        
        rms_full = np.sqrt(rms_full / max(count, 1))
        rms_free = np.sqrt(rms_free / max(count, 1))
        print(f"\nRMS error (full J_N): {rms_full:.6f}")
        print(f"RMS error (free J_N^free): {rms_free:.6f}")
        print(f"Improvement ratio: {rms_free/rms_full:.2f}x")
        
        # ========================================
        # KEY COMPARISON 2: Spectral shift vs GW prime sum
        # For test function h(E) = exp(-E^2/(2*sigma^2))
        # LHS = Tr h(J_N) - Tr h(J_N^free)
        # RHS = -int h'(E) S(E) dE  (GW prime sum)
        # ========================================
        print(f"\n--- Trace Identity Test ---")
        
        for sigma in [10.0, 30.0, 50.0]:
            def h(E):
                return np.exp(-E**2 / (2 * sigma**2))
            
            lhs = np.sum(h(eigs_full)) - np.sum(h(eigs_free))
            
            # Compute using known zeros: Tr h(J_) - Tr h(J_^free) should =  h(_k) -  h(E) '/ dE
            # For finite N, compare with partial sums
            h_zeros = np.sum(h(np.array(ZETA_ZEROS[:n_compare])))
            
            # Free trace via Weyl: Tr h(J_free) ~  h(a_k) ~  h(E) '/ dE
            h_free_trace = np.sum(h(eigs_free))
            
            print(f"  sigma={sigma:5.0f}: Tr h(J)-Tr h(J_free) = {lhs:12.8f}")
            print(f"             h(_k) -  h(_free)  = {h_zeros - h_free_trace:12.8f}")
            print(f"            Ratio: {(h_zeros - h_free_trace) / lhs:.6f}" if abs(lhs) > 1e-15 else "            Ratio: inf")
        
        # ========================================
        # KEY COMPARISON 3: Per-prime eigenvalue decomposition
        # ========================================
        print(f"\n--- Per-Prime Eigenvalue Shift (top 10 eigenvalues, first 5 primes) ---")
        a_rank1 = a_free.copy()
        a_rank1[0] += Delta
        eigs_rank1 = eigvalsh_tridiagonal(a_rank1, b_free)
        
        for p in [2, 3, 5, 7, 11]:
            dP_p, oP_p, alpha_p = build_single_prime_perturbation(N, g, p)
            a_p = a_rank1.copy()
            a_p += dP_p
            b_p = b_free + oP_p
            eigs_p = eigvalsh_tridiagonal(a_p, b_p)
            shift = eigs_p - eigs_rank1
            
            print(f"  p={p:3d} (alpha_p={alpha_p:.6f}):")
            print(f"    Shifts: [{', '.join(f'{s:8.5f}' for s in shift[:5])} ...]")
            print(f"    max|shift| = {np.max(np.abs(shift)):.6f}, "
                  f"RMS shift = {np.sqrt(np.mean(shift**2)):.6f}")
        
        # ========================================
        # KEY COMPARISON 4: Diagonal vs Off-diagonal contribution
        # This is the unexplored direction from the obstruction map
        # ========================================
        print(f"\n--- Diagonal vs Off-Diagonal Perturbation Effect ---")
        
        a_diag = a_rank1.copy()
        a_diag += dP  # diagonal only
        eigs_diag = eigvalsh_tridiagonal(a_diag, b_free)
        
        a_offdiag = a_rank1.copy()
        b_offdiag = b_free + oP  # off-diagonal only
        eigs_offdiag = eigvalsh_tridiagonal(a_offdiag, b_offdiag)
        
        err_full_vec = eigs_full[:n_compare] - np.array(ZETA_ZEROS[:n_compare])
        err_free_vec = eigs_free[:n_compare] - np.array(ZETA_ZEROS[:n_compare])
        err_diag_vec = eigs_diag[:n_compare] - np.array(ZETA_ZEROS[:n_compare])
        err_offdiag_vec = eigs_offdiag[:n_compare] - np.array(ZETA_ZEROS[:n_compare])
        
        print(f"  RMS error (free only):     {np.sqrt(np.mean(err_free_vec**2)):.6f}")
        print(f"  RMS error (diag perturb):  {np.sqrt(np.mean(err_diag_vec**2)):.6f}")
        print(f"  RMS error (offdiag pert):  {np.sqrt(np.mean(err_offdiag_vec**2)):.6f}")
        print(f"  RMS error (full perturb):  {np.sqrt(np.mean(err_full_vec**2)):.6f}")
        
        # ========================================
        # KEY COMPARISON 5: Krein determinant on energy grid
        # ========================================
        print(f"\n--- Krein Determinant Spectral Shift ---")
        E_range = np.linspace(max(17, eigs_free[0] - 1), eigs_free[-1] + 1, 100)
        xi_from_eigs = np.zeros(len(E_range))
        for i, E in enumerate(E_range):
            xi_from_eigs[i] = np.searchsorted(eigs_full, E, side='right') - np.searchsorted(eigs_free, E, side='right')
        
        eps = 0.1
        args_krein = np.zeros(len(E_range))
        for i, E in enumerate(E_range):
            z = complex(E, eps)
            Jf = np.diag(a_free - z + 0j)
            for k in range(N - 1):
                Jf[k, k+1] = b_free[k]
                Jf[k+1, k] = b_free[k]
            Pm = np.diag(dP)
            for k in range(N - 1):
                Pm[k, k+1] = oP[k]
                Pm[k+1, k] = oP[k]
            try:
                G = np.linalg.inv(Jf)
                det_val = np.linalg.det(np.eye(N) + G @ Pm)
                args_krein[i] = np.angle(det_val)
            except:
                args_krein[i] = np.nan
        
        valid = ~np.isnan(args_krein)
        if np.any(valid):
            corr = np.corrcoef(xi_from_eigs[valid], args_krein[valid] / np.pi)[0, 1]
            print(f"  Correlation between xi(eigs) and arg(det)/(pi): {corr:.6f}")
            
            # Check: does arg(det) ~ pi * (N_J(E) - N_free(E)) ?
            residual = args_krein[valid] / np.pi - xi_from_eigs[valid]
            print(f"  arg(det)/pi - xi: mean={np.mean(residual):.6f}, std={np.std(residual):.6f}")
            print(f"  (Should be ~0 for Birman-Krein to hold)")
        
        # ========================================
        # KEY COMPARISON 6: Fourier analysis of eigenvalue shifts
        # If eigenvalues track zeta zeros, the shifts _k = _k - a_k
        # should have Fourier content at frequencies {m log p}
        # ========================================
        print(f"\n--- Fourier Analysis of Eigenvalue Shifts ---")
        shifts = eigs_full - a_free  # eigenvalue - diagonal entry
        n_fft = min(N, 80)
        shifts_trunc = shifts[:n_fft]
        
        # Use the diagonal entries as "time" axis (approximately equally spaced in E)
        E_axis = a_free[:n_fft]
        
        # Compute FFT of shifts vs E
        dE = np.mean(np.diff(E_axis))
        fft_shifts = np.fft.rfft(shifts_trunc)
        freqs = np.fft.rfftfreq(n_fft, d=dE)
        
        # Convert to angular frequencies
        omega = 2 * np.pi * freqs
        
        # Find peaks
        power = np.abs(fft_shifts)
        top_indices = np.argsort(power[1:])[-10:] + 1  # skip DC
        
        print(f"  Top 10 Fourier frequencies in eigenvalue shifts:")
        print(f"  {'omega':>12} {'power':>12} {'nearest_log_p':>14} {'prime':>6}")
        for idx in top_indices:
            w = omega[idx]
            pwr = power[idx]
            # Find nearest log(p)
            best_p = 0
            best_diff = float('inf')
            for p in primerange(2, 100):
                for m in range(1, 5):
                    diff = abs(w - m * np.log(p))
                    if diff < best_diff:
                        best_diff = diff
                        best_p = p
                        best_m = m
            label = f"{best_m}*log({best_p})" if best_m > 1 else f"log({best_p})"
            print(f"  {w:12.4f} {pwr:12.6f} {np.log(best_p):14.4f} p={best_p:3d} ({label})")
    
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    print("""
The key question: do the eigenvalues of J_N (free + rank1 + Paley-Legendre)
converge to the Riemann zeta zeros as N -> infinity?

Evidence from this computation:
1. Eigenvalue-zero RMS error (with perturbation) vs (without perturbation)
2. Whether the perturbation IMPROVES the approximation
3. Whether off-diagonal terms contribute meaningfully
4. Whether the Krein determinant arg correlates with the spectral shift
5. Whether Fourier analysis of shifts shows prime-frequency peaks
""")

if __name__ == "__main__":
    main()
