# -*- coding: utf-8 -*-
"""Table 1 of the paper (bump data), computed at 60 digits.

Prints, for every bump of F4, F5, F6: the peak u_* = p/(p+q), the amplitude
B(u_*), and the value at 0.5761; plus the mark-perturbation data used in the
pointwise certification (largest sum_j B_j(e_i) vs |F1(e_i)|).

Run:  python scripts/certify_table.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mpmath import mp, mpf
from fewnomial_ladder.functions import F1, bump, BUMP_DEFS, MARKS

mp.dps = 60

print("=" * 72)
print("Table 1: bump data (60 digits)")
print("=" * 72)
print("instance | peak u_*        | amplitude B(u_*) | value at 0.5761")
print("-" * 72)
for t in (4, 5, 6):
    for amp, p, q in BUMP_DEFS[t]:
        u_star = mpf(p) / (p + q)
        print("F%d       | %s | %s | %s"
              % (t, mp.nstr(u_star, 8), mp.nstr(bump(amp, p, q, u_star), 8),
                 mp.nstr(bump(amp, p, q, mpf("0.5761")), 8)))

print()
print("mark perturbations (pointwise certification of Lemma 'multi'):")
print("mark e_i | F1(e_i)         | max_t sum_j B_j(e_i) | sign preserved?")
print("-" * 72)
for e in MARKS:
    v = F1(e)
    pert = {}
    for t in (4, 5, 6):
        pert[t] = sum(bump(amp, p, q, e) for amp, p, q in BUMP_DEFS[t])
    worst = max(pert.values())
    print("%s | %s | %s (F%d) | %s"
          % (mp.nstr(e, 8), mp.nstr(v, 10), mp.nstr(worst, 8),
             max(pert, key=pert.get),
             "YES" if worst < abs(v) else "NO"))
    assert worst < abs(v)

print()
print("K-margin mu = -F1(0.55) =", mp.nstr(-F1(mpf("0.55")), 10),
      " (F1 increasing on (0, 1/2))")
print("cluster margin m = 2.6e-4;  |F1(0.55)| =", mp.nstr(abs(F1(mpf("0.55"))), 10))
