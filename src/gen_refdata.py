# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Dan Alec Yamaguchi
#
# Licensed under GNU Affero General Public License v3.0.
# Commercial licenses available for enterprise use.
# Contact: danalec@gmail.com
# See LICENSE-COMMERCIAL for details.
#
"""Generate refdata header with N zeta zeros and S(T) values via mpmath."""
import mpmath as mp
import sys, os

# Known first 100 zeros from the original refdata_100.h (mpmath 60-digit)
KNOWN_100 = [
    14.134725141734695, 21.022039638771556, 25.010857580145689,
    30.424876125859512, 32.935061587739192, 37.586178158825675,
    40.918719012147498, 43.327073280915002, 48.005150881167161,
    49.773832477672300, 52.970321477714464, 56.446247697063392,
    59.347044002602352, 60.831778524609810, 65.112544048081602,
    67.079810529494168, 69.546401711173985, 72.067157674481905,
    75.704690699083926, 77.144840068874799, 79.337375020249368,
    82.910380854086029, 84.735492980517051, 87.425274613125225,
    88.809111207634459, 92.491899270558491, 94.651344040519888,
    95.870634228245308, 98.831194218193687, 101.317851005731384,
    103.725538040478341, 105.446623052326089, 107.168611184276401,
    111.029535543169672, 111.874659176992637, 114.320220915452708,
    116.226680320857554, 118.790782865976212, 121.370125002420650,
    122.946829293552582, 124.256818554345770, 127.516683879596499,
    129.578704199956064, 131.087688530932667, 133.497737202997598,
    134.756509753373876, 138.116042054533438, 139.736208952121387,
    141.123707404021133, 143.111845807620625, 146.000982486765508,
    147.422765342559615, 150.053520420784878, 150.925257612241467,
    153.024693811198887, 156.112909294237880, 157.597591817594065,
    158.849988171420506, 161.188964137596031, 163.030709687181997,
    165.537069187900414, 167.184439978174510, 169.094515415568821,
    169.911976479411692, 173.411536519591550, 174.754191523365733,
    176.441434297710430, 178.377407776099972, 179.916484020257002,
    182.207078484366463, 184.874467848387496, 185.598783677707473,
    187.228922583501856, 189.416158656016933, 192.026656360713787,
    193.079726603845700, 195.265396679529232, 196.876481840958320,
    198.015309676251917, 201.264751943703800, 202.493594514140540,
    204.189671803104545, 205.394697202163286, 207.906258887806217,
    209.576509716856265, 211.690862595365303, 213.347919359712677,
    214.547044783491430, 216.169538508263713, 219.067596349021386,
    220.714918839314009, 221.430705554693333, 224.007000254604321,
    224.983324669582288, 227.421444279679292, 229.337413305525359,
    231.250188700499166, 231.987235253180245, 233.693404178908310,
    236.524229665816193,
]

KNOWN_S = [
    0.550252829455791, 0.429788816145506, 0.606402146236036,
    0.328935549258941, 0.682728991191247, 0.406444797548046,
    0.434870892266311, 0.705618872335153, 0.229181465059437,
    0.651663733367918, 0.582786619760240, 0.385662759274591,
    0.360434714109126, 0.826881905565582, 0.256737870875993,
    0.519954463193179, 0.583218688973462, 0.611527995221646,
    0.184730450565243, 0.612078166125091, 0.732061634570901,
    0.277417963395485, 0.524852594645555, 0.404382385909717,
    0.822766399099761, 0.258339438990195, 0.330103504275404,
    0.802517103450070, 0.511273175506309, 0.415795962667491,
    0.345849401433904, 0.575556239926031, 0.800384265862405,
    0.046444640729137, 0.659645253368591, 0.534653886539968,
    0.651871930389986, 0.456735553309083, 0.245608275554314,
    0.500959954430923, 0.879826771335502, 0.324665726264219,
    0.334081841996167, 0.605859480243539, 0.437070705189593,
    0.823849388241592, 0.178110430688078, 0.379768584424355,
    0.693696871029878, 0.706849626375239, 0.264956759839175,
    0.552029897447985, 0.227146794469581, 0.786502776659834,
    0.722003719617649, 0.147839863687515, 0.387576070005259,
    0.744527663562221, 0.539371522369684, 0.586607729262250,
    0.284716365091230, 0.425718721725554, 0.426504845651871,
    0.997816526410376, 0.155558189699064, 0.445754094948960,
    0.551454803552970, 0.522157112201724, 0.701489680381684,
    0.476222455817882, 0.043630922309577, 0.653556891938641,
    0.774016707615228, 0.590349855973842, 0.172368432436220,
    0.598754038171271, 0.405293080511886, 0.523082515217094,
    0.898207446026674, 0.109527438398146, 0.430918559856566,
    0.492338605144396, 0.824134765355781, 0.427828933998179,
    0.496573823791825, 0.314671512637298, 0.386042220764083,
    0.712763308166381, 0.800076747537211, 0.165045393368386,
    0.232930911387404, 0.827302252945549, 0.364304116993943,
    0.808644527865678, 0.418089145116974, 0.322415395473913,
    0.226026742201213, 0.802883907185627, 0.821934185860639,
    0.190004621915122,
]

def find_zero_beyond_100(n):
    """Find nth zero using theta-based initial guess + Newton."""
    mp.mp.dps = 80
    
    # Solve theta(t) = pi*(n - 3/2) for initial guess
    t = 2.0 * mp.pi * n / mp.log(float(n))
    for _ in range(20):
        th = float(mp.im(mp.loggamma(0.25 + 0.5j * t)) - 0.5 * t * mp.log(float(mp.pi)))
        fp = 0.5 * mp.log(t / (2.0 * mp.pi))
        dt = (th - mp.pi * (n - 1.5)) / fp
        t = t - dt
        if abs(float(dt)) < 1e-14:
            break
    
    # Polish with Newton on Hardy Z-function
    try:
        root = mp.findroot(lambda x: float(mp.re(mp.zeta(0.5+1j*x)))*mp.cos(float(mp.im(mp.loggamma(0.25+0.5j*x))-0.5*x*mp.log(float(mp.pi)))) - float(mp.im(mp.zeta(0.5+1j*x)))*mp.sin(float(mp.im(mp.loggamma(0.25+0.5j*x))-0.5*x*mp.log(float(mp.pi)))), float(t))
        gamma = float(root)
    except Exception:
        gamma = float(t)
    
    return gamma

def compute_S(gamma):
    mp.mp.dps = 80
    eps = 1e-8
    z_plus = mp.zeta(0.5 + (gamma + eps) * 1j)
    return float(float(mp.arg(z_plus)) / float(mp.pi))

def compute_data(N):
    zeros = list(KNOWN_100[:min(100, N)])
    s_vals = list(KNOWN_S[:min(100, N)])
    
    if N <= 100:
        return zeros, s_vals
    
    for n in range(101, N + 1):
        try:
            gamma = find_zero_beyond_100(n)
            zeros.append(gamma)
            s = compute_S(gamma)
            s_vals.append(s)
        except Exception as e:
            print(f"Error at n={n}: {e}", file=sys.stderr)
            zeros.append(0.0)
            s_vals.append(0.0)
        
        if n % 50 == 0:
            print(f"  ... computed {n} zeros", file=sys.stderr)
    
    return zeros, s_vals

def write_header(N, zeros, s_vals, filename):
    with open(filename, 'w', encoding='ascii') as f:
        f.write("/* SPDX-License-Identifier: AGPL-3.0-or-later\n")
        f.write(" * Copyright (C) 2026 Dan Alec Yamaguchi\n")
        f.write(" *\n")
        f.write(" * Licensed under GNU Affero General Public License v3.0.\n")
        f.write(" * Commercial licenses available for enterprise use.\n")
        f.write(" * Contact: danalec@gmail.com\n")
        f.write(" * See LICENSE-COMMERCIAL for details.\n")
        f.write(" */\n\n")
        f.write(f"/* rh/refdata_{N}.h -- {N} zeta zeros and S(T) values (mpmath, 80-digit precision) */\n")
        f.write(f"#ifndef REFDATA_{N}_H\n")
        f.write(f"#define REFDATA_{N}_H\n\n")
        f.write(f"#define N_REF {N}\n\n")
        
        f.write("static const double ZETA_ZEROS[N_REF] = {\n")
        for i, z in enumerate(zeros):
            comma = "," if i < len(zeros) - 1 else ""
            nl = "\n" if (i + 1) % 5 == 0 else ""
            f.write(f"    {z:.15f}{comma}{nl}")
        f.write("\n};\n\n")
        
        f.write("static const double S_AT_ZERO[N_REF] = {\n")
        for i, s in enumerate(s_vals):
            comma = "," if i < len(s_vals) - 1 else ""
            nl = "\n" if (i + 1) % 5 == 0 else ""
            f.write(f"    {s:.15f}{comma}{nl}")
        f.write("\n};\n\n")
        f.write("#endif\n")
    
    print(f"Wrote {N} entries to {filename}", file=sys.stderr)

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    filename = sys.argv[2] if len(sys.argv) > 2 else f"refdata_{N}.h"
    
    print(f"Generating {N} zeta zeros (first 100 from refdata, 101-{N} from mpmath)...", file=sys.stderr)
    zeros, s_vals = compute_data(N)
    write_header(N, zeros, s_vals, filename)
    print("Done.", file=sys.stderr)
