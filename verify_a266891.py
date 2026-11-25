"""Verification for "Weighted partitions with bounded multiplicities" (OEIS A266891).

a(n) = [x^n] prod_{k>=1} (1 + k x^k)^k.

1. Exact coefficients a(0..NMAX) by integer arithmetic; compared with the OEIS
   data and, if available, with the b-file b266891.txt (downloaded unless
   --no-download is given).
2. Table: log a(n), the explicit main term, log M_K, the saddle value
   min_s Phi_n(s), and the Gaussian-corrected value min Phi - (1/2) log(2 pi V).

Usage: python3 verify_a266891.py [NMAX] [--no-download]
"""
import os
import sys
import urllib.request
from math import comb, exp, log, log1p, pi, sqrt

OEIS = [1, 1, 4, 13, 29, 81, 188, 456, 1030, 2405, 5295, 11611, 25246, 53552, 113332,
        235685, 486011, 990840, 2006567, 4018010, 7992003, 15768511, 30875424,
        60060509, 116042548, 222817961, 425200270, 806991037, 1522748592, 2858792520,
        5339457208, 9924370365]

args = [a for a in sys.argv[1:] if not a.startswith('--')]
NMAX = int(args[0]) if args else 2000

# ------------------------------------------------------------ 1. coefficients
a = [1] + [0] * NMAX
for k in range(1, NMAX + 1):
    terms = [(j * k, comb(k, j) * k**j) for j in range(1, min(k, NMAX // k) + 1)]
    for i in range(NMAX, -1, -1):
        s = a[i]
        for d, c in terms:
            if d > i:
                break
            s += c * a[i - d]
        a[i] = s
assert a[:len(OEIS)] == OEIS
print("1. coefficients agree with the OEIS data (n=0..%d)" % (len(OEIS) - 1))

if not os.path.exists('b266891.txt') and '--no-download' not in sys.argv:
    try:
        req = urllib.request.Request('https://oeis.org/A266891/b266891.txt',
                                     headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=60) as f, open('b266891.txt', 'wb') as g:
            g.write(f.read())
    except Exception as exc:                     # network is optional
        print("   (b-file download failed: %s)" % exc)
if os.path.exists('b266891.txt'):
    checked = 0
    for line in open('b266891.txt'):
        line = line.strip()
        if line and not line.startswith('#'):
            i, v = map(int, line.split()[:2])
            if i <= NMAX:
                assert a[i] == int(v), i
                checked += 1
    print("   coefficients agree with the b-file for %d indices" % checked)


# ------------------------------------------------------------ 2. table
def ln(x):
    b = x.bit_length()
    sh = max(0, b - 60)
    return log(x >> sh) + sh * log(2)


KMAX = 50000


def phi_mv(s, n):
    P, m, v = s * n, 0.0, 0.0
    for k in range(1, KMAX):
        x = k * exp(-s * k)
        if x < 1e-300:
            break
        p = x / (1 + x)
        P += k * log1p(x)
        m += k * k * p
        v += k**3 * p * (1 - p)
    return P, m, v


print("2.    n   K   log a(n)  main term   log M_K   min Phi  min Phi-0.5log(2piV)  diff")
for n in (100, 200, 500, 1000, 1500, 2000, 3000, 5000, 10000):
    if n > NMAX:
        break
    lo, hi = 1e-5, 3.0
    for _ in range(100):
        s = (lo + hi) / 2
        _, m, _ = phi_mv(s, n)
        lo, hi = (s, hi) if m > n else (lo, s)
    P, m, v = phi_mv(s, n)
    K = 0
    while (K + 1) * (K + 2) * (2 * K + 3) // 6 <= n:
        K += 1
    lM = sum(k * log(k) for k in range(1, K + 1))
    main = n**(2 / 3) * (2 * log(3 * n) - 3) / (4 * 3**(1 / 3))
    g = P - 0.5 * log(2 * pi * v)
    la = ln(a[n])
    print("   %5d %3d %9.3f %9.3f %9.3f %9.3f %12.4f %12.4f" % (n, K, la, main, lM, P, g, la - g))
print("ALL CHECKS PASSED")
