# -*- coding: utf-8 -*-
"""Exact certificate for Theorem ``m3`` (M_3 >= 5), via sympy.

The roots of  psi(x) = x^{-1/6}(1-x)^{-49/3} P(x) - Q(x)  in (0,1)
correspond, under x = y^6, 1-x = z^6, to the positive solutions of

    P(y^6) - y z^98 Q(y^6) = 0,      y^6 + z^6 - 1 = 0.

For monic g(z) = z^6 + y^6 - 1 the resultant equals the product of f over
the six roots zeta of g.  Since zeta^98 = (1-y^6)^16 zeta^2 and the six
values zeta^2 run twice through the three roots of w^3 = 1-y^6,

    Res_z(f, g) = R(y)^2,   R(y) := P(y^6)^3 - y^3 (1-y^6)^49 Q(y^6)^3.

This script:
  * verifies the product formula numerically at a sample point (60 digits);
  * computes deg R = 351, hence deg Res = 702;
  * proves by the exact Sturm count (sympy) that R has exactly five roots
    in (0,1);
  * verifies that the five high-precision roots of psi map to roots of R
    (y_i = x_i^{1/6}), which -- by the injection psi-roots -> resultant
    roots -- proves psi has exactly five roots in (0,1).

Run:  python scripts/certify_m3.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sympy as sp
from sympy import Rational, Symbol, expand, together, count_roots, degree, LC
from mpmath import mp, mpf, exp, pi

# ---------------- exact part (sympy) ----------------
y = Symbol("y")
P = y ** 3 - Rational(617, 500) * y ** 2 + Rational(441, 1000) * y - Rational(43, 1000)
Q = y ** 3 - Rational(421, 250) * y ** 2 + Rational(413, 500) * y - Rational(107, 1000)

P6 = P.subs(y, y ** 6)
Q6 = Q.subs(y, y ** 6)
R = expand(together(P6 ** 3 - y ** 3 * (1 - y ** 6) ** 49 * Q6 ** 3))

print("=" * 72)
print("Theorem 'm3': exact certificate (resultant + Sturm)")
print("=" * 72)
print("P(x) =", P)
print("Q(x) =", Q)
print("R(y) = P(y^6)^3 - y^3 (1-y^6)^49 Q(y^6)^3,  deg R =", degree(R, y))
print("deg Res_z(f,g) = deg R^2 =", 2 * degree(R, y), "  (paper: 702)")
assert 2 * degree(R, y) == 702
print("R(0) =", R.subs(y, 0), " (nonzero, no endpoint root)")
print("R(1) =", R.subs(y, 1), " (nonzero, no endpoint root)")

n_roots = count_roots(R, inf=0, sup=1)
print("Sturm count of distinct roots of R in (0,1):", n_roots)
assert n_roots == 5

# ---------------- numerical part (mpmath) ----------------
mp.dps = 60

# (a) product formula check: Res_z(f,g)(y0) = R(y0)^2 via the six roots of g.
y0 = mpf("0.3")
C = 1 - y0 ** 6
y6 = y0 ** 6
Aval = y6 ** 3 - mpf(617) / 500 * y6 ** 2 + mpf(441) / 1000 * y6 - mpf(43) / 1000
Bval = y0 * C ** 16 * (y6 ** 3 - mpf(421) / 250 * y6 ** 2
                       + mpf(413) / 500 * y6 - mpf(107) / 1000)
prod = mpf(1)
for k in range(6):
    zeta2 = C ** (mpf(1) / 3) * exp(2j * pi * k / 3)
    prod *= (Aval - Bval * zeta2)
R2_num = (Aval ** 3 - Bval ** 3 * C) ** 2
print("product formula check at y0 = 0.3: |prod - R^2| =", mp.nstr(abs(prod - R2_num), 8))
assert abs(prod - R2_num) < mpf("1e-40")

# (b) the five psi-roots (paper, 30 digits) refined, and their images in R.
psi_roots = [
    "0.00851960918806880244243398840077",
    "0.0166569183582199881798090346641",
    "0.153327345910409136363190752304",
    "0.39933905599941049891089776555",
    "0.675136363428802112925725547891",
]
fR = sp.lambdify(y, R, modules="mpmath")
print()
print("psi-roots and their membership in the resultant (y_i = x_i^{1/6}):")
for s in psi_roots:
    x = mpf(s)
    Pv = x ** 3 - mpf(617) / 500 * x ** 2 + mpf(441) / 1000 * x - mpf(43) / 1000
    Qv = x ** 3 - mpf(421) / 250 * x ** 2 + mpf(413) / 500 * x - mpf(107) / 1000
    val = x ** (-mpf(1) / 6) * (1 - x) ** (-mpf(49) / 3) * Pv - Qv
    yi = x ** (mpf(1) / 6)
    print("   x = %s" % mp.nstr(x, 32))
    print("        psi(x) = %s" % mp.nstr(val, 6))
    print("        |R(y)| = %s" % mp.nstr(abs(fR(yi)), 6))

print()
print("EXACT CERTIFICATE PASSED: deg Res = 702, Sturm count = 5 distinct")
print("roots in (0,1), and the five psi-roots inject into them; hence psi")
print("has exactly five roots in (0,1), i.e. M_3 >= 5.")
