# -*- coding: utf-8 -*-
"""Symbolic check of Lemma ``reduction`` applied to the Mueller--Regensburger
system (eq. (mr)):

    x^5/y + b y - 1 = 0,     y^5/x + b x - 1 = 0.

With u = b x and v = y^5/x the second equation becomes u + v = 1 (case (i)
of the lemma), and substituting v = 1 - u into the first equation gives
exactly F1(u) = b^{-24/5} u^{24/5} (1-u)^{-1/5}
                + b^{4/5} u^{1/5} (1-u)^{1/5} - 1.

Run:  python scripts/verify_reduction.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sympy as sp

b, u, v = sp.symbols("b u v", positive=True)
x = u / b
y = (v * u / b) ** sp.Rational(1, 5)

eq1 = x ** 5 / y + b * y - 1        # first MR equation in (u, v)
eq2 = v + u - 1                     # second MR equation in (u, v)

F1_expr = (b ** (-sp.Rational(24, 5)) * u ** (sp.Rational(24, 5))
           * (1 - u) ** (-sp.Rational(1, 5))
           + b ** (sp.Rational(4, 5)) * u ** (sp.Rational(1, 5))
           * (1 - u) ** (sp.Rational(1, 5)) - 1)

diff = sp.simplify(eq1.subs(v, 1 - u) - F1_expr)

print("substituting v = 1 - u into the first MR equation gives F1(u):")
print("   difference =", diff)
assert diff == 0
print("CHECK PASSED: eq. (mr) reduces exactly to F1(u) = 0 with u = b x,"
      " v = y^5/x.")
