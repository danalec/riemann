# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dan Alec Yamaguchi
#
# Licensed under GNU Affero General Public License v3.0.
# Commercial licenses available for enterprise use.
# Contact: danalec@gmail.com
# See LICENSE-COMMERCIAL for details.
#
"""
Check pointwise convergence: does lambda_k^(N) -> gamma_k as N -> infinity
for FIXED k? This is the paper's key claim (line 1063).
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

ZETA_ZEROS = [
    14.134725141734693, 21.022039638771554, 25.010857580145688,
    30.424876125859513, 32.935061587739190, 37.586178158825671,
    40.918719012147495, 43.327073280914999, 48.005150881167159,
    49.773832477672302, 52.970321477714460, 56.446247697063395,
    59.347044002602903, 60.831778524609809, 65.112544048081607,
    67.079810529501173, 69.546401711173979, 72.103742237830137,
    75.704690699083933, 77.144840068874805,
]

def compute_eigenvalues(N):
    g = np.array([gram_point(n) for n in range(N + 1)])
    Lambda = g[N - 1]
    
    a = np.zeros(N)
    for n in range(N):
        gn = g[0] if n == 0 else g[n - 1]
        a[n] = gn + np.pi / np.log(gn / (2 * np.pi))
    
    b = np.zeros(N - 1)
    for n in range(N - 1):
        b[n] = np.sqrt(g[n + 1] - g[n]) * theta_prime(g[n + 1])
    
    dP = np.zeros(N)
    oP = np.zeros(N - 1)
    for p in primerange(2, int(Lambda) + 1):
        alpha_p = -np.log(p) / (2 * np.pi * np.sqrt(p))
        logp = np.log(p)
        for n in range(N):
            dP[n] += alpha_p * np.sin(logp * g[n])
        for n in range(N - 1):
            oP[n] += alpha_p * np.cos(logp * (g[n] + g[n + 1]) / 2)
    
    gamma_1 = ZETA_ZEROS[0]
    Delta = gamma_1 - (g[0] + np.pi / np.log(g[0] / (2 * np.pi)))
    
    a[0] += Delta
    a += dP
    b += oP
    
    return eigvalsh_tridiagonal(a, b)

print("Pointwise convergence check: lambda_k^(N) for fixed k as N grows")
print("=" * 80)

Ns = [30, 40, 50, 60, 80, 100, 150, 200]
k_check = [1, 2, 3, 5, 8, 10, 12, 15, 19]

# Store results
results = {}
for N in Ns:
    print(f"Computing N={N}...", end=" ", flush=True)
    eigs = compute_eigenvalues(N)
    for k in k_check:
        if k < len(eigs) and k < len(ZETA_ZEROS):
            key = k
            if key not in results:
                results[key] = []
            results[key].append((N, eigs[k] - ZETA_ZEROS[k]))
    print("done")

print()
print(f"{'k':>4} {'gamma_k':>12}", end="")
for N in Ns:
    print(f"  N={N:>4}", end="")
print()
print("-" * 120)

for k in k_check:
    if k in results:
        print(f"{k:4d} {ZETA_ZEROS[k]:12.4f}", end="")
        for N in Ns:
            val = None
            for n, err in results[k]:
                if n == N:
                    val = err
                    break
            if val is not None:
                print(f"  {val:8.4f}", end="")
            else:
                print(f"  {'---':>8}", end="")
        print()

print()
print("Trend analysis: Is the error decreasing for each fixed k?")
for k in k_check:
    if k in results and len(results[k]) >= 3:
        errs = [err for _, err in results[k]]
        Ns_k = [n for n, _ in results[k]]
        if len(errs) >= 4:
            # Fit linear trend
            coeffs = np.polyfit(Ns_k, errs, 1)
            slope = coeffs[0]
            print(f"  k={k:3d}: slope = {slope:+.6f} per unit N ({'CONVERGING' if slope * Ns_k[-1] < 0 and abs(slope) > 1e-5 else 'NOT CONVERGING'})")
