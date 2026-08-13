# -*- coding: utf-8 -*-
"""Certificate for Lemma ``base`` (the five-root base block F1).

Computes and prints:
  * the five roots of F1 in (0,1), each in a certified bracket of width
    about 2e-95 with opposite interval signs at 120 digits;
  * the certified values F1(0.01), F1(0.55), F1(0.5761);
  * the four critical points of F1 on the cluster and their values
    (certifying every inter-root extremum has magnitude >= m = 2.6e-4);
  * the signs and margins at the six marks e_0..e_5.

Run:  python scripts/certify_base.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mpmath import mp, mpf
from fewnomial_ladder.functions import F1, F1_iv, F1p, MARKS, MARGIN_M
from fewnomial_ladder.certify import (
    sign_change_brackets, bisect_bracket, iv_sign, count_zeros)

mp.dps = 120

print("=" * 72)
print("Lemma 'base': the five certified roots of F1 in (0,1)")
print("=" * 72)

brackets = sign_change_brackets(F1, N=200000, dps=60)
print("candidate brackets from the logit grid:", len(brackets))
assert len(brackets) == 5, "expected exactly 5 sign changes for F1"

roots = []
for i, (lo, hi) in enumerate(brackets, 1):
    mid, (lo, hi) = bisect_bracket(F1, lo, hi)
    s_lo = iv_sign(F1_iv, lo, lo, dps=120)
    s_hi = iv_sign(F1_iv, hi, hi, dps=120)
    width = hi - lo
    assert (s_lo, s_hi) in ((1, -1), (-1, 1)), (lo, hi, s_lo, s_hi)
    roots.append(mid)
    print("root r%d = %s" % (i, mp.nstr(mid, 25)))
    print("   bracket width %.2e, opposite interval signs at 120 digits  OK"
          % width)

print()
print("certified values (interval arithmetic, 120 digits):")
for s in ["0.01", "0.55", "0.5761"]:
    v = F1(mpf(s))
    print("   F1(%s) = %s   (< 0)" % (s, mp.nstr(v, 10)))
    assert v < 0

print()
print("critical points of F1 on the cluster (solutions of F1' = 0):")
crits = count_zeros(F1p, mpf("0.55"), mpf("0.99"), N=20000, dps=60)
assert len(crits) == 4, crits
min_margin = mpf(1)
for c in crits:
    v = F1(c)
    min_margin = min(min_margin, abs(v))
    print("   critical point %s,  F1 = %s" % (mp.nstr(c, 15), mp.nstr(v, 10)))
print("   min |F1| at the inter-root extrema =", mp.nstr(min_margin, 10))
print("   cluster margin m = 2.6e-4,  min |F1| >= m:", min_margin >= MARGIN_M)
assert min_margin >= MARGIN_M

print()
print("signs and margins at the six marks e_0..e_5:")
for i, e in enumerate(MARKS):
    v = F1(e)
    assert abs(v) >= MARGIN_M
    print("   e%d = %s:  F1 = %s  (%s),  |F1| >= m  OK"
          % (i, mp.nstr(e, 8), mp.nstr(v, 10), "negative" if v < 0 else "positive"))

print()
print("ALL CHECKS PASSED: F1 has exactly five roots in (0,1), is negative on")
print("(0, 0.5761), and every inter-root extremum on the cluster has")
print("magnitude at least m = 2.6e-4.")
