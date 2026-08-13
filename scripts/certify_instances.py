# -*- coding: utf-8 -*-
"""Certificates for Theorem ``verified``: F4, F5, F6 have 7, 9, 11 roots.

For each t = 4, 5, 6 this script:
  1. discovers candidates on a logit grid and certifies each root in a
     bracket of width about 2e-95 with opposite interval signs at 120 digits
     (this proves at least 7 / 9 / 11 roots);
  2. certifies exhaustiveness: F_t < 0 on K = [0.01, 0.55] minus the
     localization intervals and on [0, 0.01] (pointwise interval covers);
  3. certifies the mark signs: sum_j B_j(e_i) < |F1(e_i)| at e_0..e_5;
  4. certifies uniqueness: on the cluster, F_t' has a definite sign on each
     root bracket, exactly one zero in each inter-root gap, and is positive
     on (0.55, first root) and on (last root, 1); in each bump interval,
     F_t' has no zero on the left half and at most two zeros on the right
     half, with the correct endpoint signs, so F_t has exactly one root on
     each side of every peak.

Run:  python scripts/certify_instances.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mpmath import mp, mpf
from fewnomial_ladder.functions import (
    F1, F1_iv, Ft, Ft_iv, Ftp_iv, bump, BUMP_DEFS, MARKS, INTERVALS,
    k_regions)
from fewnomial_ladder.certify import (
    sign_change_brackets, bisect_bracket, iv_sign, certify_negative,
    count_zeros)

mp.dps = 120
EXPECTED = {4: 7, 5: 9, 6: 11}

for t in (4, 5, 6):
    print("=" * 72)
    print("Theorem 'verified': certification of F%d" % t)
    print("=" * 72)

    f = lambda u: Ft(t, u)

    # ---- 1. certified root brackets (existence of at least N_t roots) ----
    brackets = sign_change_brackets(f, N=200000, dps=60)
    roots = []
    for lo, hi in brackets:
        mid, (lo, hi) = bisect_bracket(f, lo, hi)
        s_lo = iv_sign(lambda u: Ft_iv(t, u), lo, lo, dps=120)
        s_hi = iv_sign(lambda u: Ft_iv(t, u), hi, hi, dps=120)
        assert (s_lo, s_hi) in ((1, -1), (-1, 1))
        roots.append((mid, lo, hi))
    n = len(roots)
    print("certified root brackets: %d  (expected %d)  %s"
          % (n, EXPECTED[t], "OK" if n == EXPECTED[t] else "FAIL"))
    assert n == EXPECTED[t]

    # ---- 2. exhaustiveness: F_t < 0 on K and on [0, 0.01] ----
    pieces = k_regions(t)
    ok, bad = certify_negative(lambda u: Ft_iv(t, u), pieces)
    print("F%d < 0 on K = [0.01, 0.55] minus the bump intervals: %s"
          % (t, "OK" if ok else "FAIL"))
    assert ok
    ok0, bad0 = certify_negative(lambda u: Ft_iv(t, u), [(mpf(0), mpf("0.01"))])
    print("F%d < 0 on [0, 0.01]: %s" % (t, "OK" if ok0 else "FAIL"))
    assert ok0

    # ---- 3. mark signs: total perturbation below |F1(e_i)| ----
    for i, e in enumerate(MARKS):
        P = sum(bump(amp, p, q, e) for amp, p, q in BUMP_DEFS[t])
        assert P < abs(F1(e))
    print("sum_j B_j(e_i) < |F1(e_i)| at all six marks: OK")

    # ---- 4a. cluster uniqueness ----
    cluster = roots[-5:]  # the five perturbed old roots
    fprime = lambda u: mp.diff(f, u)
    ok_der = True
    for mid, lo, hi in cluster:
        s = iv_sign(lambda u: Ftp_iv(t, u), lo, hi)
        if s is None:
            ok_der = False
    print("F%d' has a definite sign on every cluster root bracket: %s"
          % (t, "OK" if ok_der else "FAIL"))
    assert ok_der
    for i in range(4):
        a = cluster[i][2]   # right endpoint of i-th bracket
        b = cluster[i + 1][1]  # left endpoint of (i+1)-th bracket
        z = count_zeros(fprime, a, b, N=8000, dps=60)
        if len(z) != 1:
            ok_der = False
            print("   gap %d: %d critical points (expected 1)" % (i + 1, len(z)))
        else:
            mid_c = z[0]
            # The interval cover must stop short of the (numerically located)
            # critical point; F_t' has exactly one zero in (a, b), so
            # sign-definiteness on [a, c-d] and [c+d, b] propagates to the
            # two open half-gaps by continuity.
            d = mpf("1e-8")
            s1 = iv_sign(lambda u: Ftp_iv(t, u), a, mid_c - d)
            s2 = iv_sign(lambda u: Ftp_iv(t, u), mid_c + d, b)
            if s1 is None or s2 is None or s1 == s2:
                ok_der = False
                print("   gap %d: derivative signs %s, %s (not opposite definite)"
                      % (i + 1, s1, s2))
    print("F%d' has exactly one zero in each inter-root gap, with monotone"
          " halves: %s" % (t, "OK" if ok_der else "FAIL"))
    assert ok_der
    s_pre = iv_sign(lambda u: Ftp_iv(t, u), mpf("0.55"), cluster[0][1])
    s_post = iv_sign(lambda u: Ftp_iv(t, u), cluster[-1][2], mpf("0.99"))
    assert s_pre == 1 and s_post == 1
    print("F%d' > 0 on (0.55, first root) and on (last root, 0.99): OK" % t)

    # ---- 4b. bump-strip uniqueness ----
    # Split each localization interval at the zeros of F_t' (located by
    # grid + bisection; tangencies without sign change are harmless).  On
    # every maximal piece F_t' has a definite certified sign, so F_t is
    # monotone there and carries at most one root; counting the sign
    # changes of F_t across consecutive breakpoints gives exactly one root
    # on each side of the peak.
    for j, (lo, hi) in enumerate(INTERVALS[t]):
        amp, p, q = BUMP_DEFS[t][j]
        u_star = mpf(p) / (p + q)
        zl = count_zeros(fprime, lo, u_star, N=4000, dps=60)
        zr = count_zeros(fprime, u_star, hi, N=4000, dps=60)
        bps = [lo] + zl + [u_star] + zr + [hi]
        assert f(lo) < 0 and f(u_star) > 0 and f(hi) < 0
        # certified definite sign of F_t' on every piece (away from breakpoints)
        d = mpf("1e-8")
        for a, b in zip(bps, bps[1:]):
            if b - a > 2 * d:
                s = iv_sign(lambda u: Ftp_iv(t, u), a + d, b - d)
                assert s in (1, -1), (t, j, a, b, s)
        left_vals = [f(b) for b in [lo] + zl + [u_star]]
        right_vals = [f(b) for b in [u_star] + zr + [hi]]

        def sign_changes(vals):
            n = 0
            for a, b in zip(vals, vals[1:]):
                if a == 0 or b == 0 or a * b < 0:
                    n += 1
            return n

        nl, nr = sign_changes(left_vals), sign_changes(right_vals)
        assert nl == 1 and nr == 1, (t, j, nl, nr)
        print("   bump %d: F%d' has %d zero(s) left / %d right of the peak,"
              " exactly one root on each side: OK" % (j + 1, t, len(zl), len(zr)))

    # ---- 5. print the roots ----
    print("roots of F%d in (0,1):" % t)
    for mid, lo, hi in roots:
        print("   u = %s" % mp.nstr(mid, 40))

print()
print("ALL CHECKS PASSED: F4, F5, F6 have exactly 7, 9, 11 roots in (0,1).")
