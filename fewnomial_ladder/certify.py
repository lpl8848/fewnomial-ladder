# -*- coding: utf-8 -*-
"""Shared certification machinery (grid discovery, bisection, interval signs).

Implements the protocol of Section "Certification protocol" of the paper:

1. discovery:  logit grid u = 1/(1+e^{-t}), t in [-40,40], sign changes
   are candidate root locations (candidate generation only);
2. certification: each candidate is refined by bisection to a bracket of
   width about 2e-95 around a root, and the sign change across the bracket
   is re-verified with mpmath interval arithmetic at 120 digits;
3. exhaustiveness/uniqueness: definite signs are certified on interval
   covers (mpmath midpoint-radius intervals; see README for the exactness
   caveat -- the rational-data count of Theorem m3 is certified exactly
   with sympy instead).
"""

from mpmath import mp, iv, mpf, exp


def logit_grid(N, t_lo=-40, t_hi=40):
    """N+1 logit-spaced points in (0,1)."""
    return [1 / (1 + exp(-(t_lo + (t_hi - t_lo) * i / N))) for i in range(N + 1)]


def sign_change_brackets(f, N=200000, dps=60):
    """Scan the logit grid for sign changes of f; returns [(lo, hi), ...]."""
    mp.dps = dps
    brackets = []
    prev_u = prev_v = None
    for u in logit_grid(N):
        v = f(u)
        if prev_v is not None and ((v < 0) != (prev_v < 0)):
            brackets.append((prev_u, u))
        prev_u, prev_v = u, v
    return brackets


def bisect_bracket(f, lo, hi, width=mpf("2e-95"), dps=120, maxiter=4000):
    """Bisect a sign-change bracket to width <= 'width' at dps digits.

    Requires f(lo) * f(hi) < 0.  Returns (mid, (lo, hi)).
    """
    mp.dps = dps
    lo, hi = mpf(lo), mpf(hi)
    flo, fhi = f(lo), f(hi)
    if not (flo * fhi < 0):
        raise ValueError("bracket does not carry a sign change")
    for _ in range(maxiter):
        mid = (lo + hi) / 2
        if hi - lo <= width:
            return mid, (lo, hi)
        fm = f(mid)
        if fm == 0:
            return mid, (mid, mid)
        if flo * fm < 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    raise RuntimeError("bisection did not converge")


def _iv_interval(lo, hi):
    """Build the mpi [lo, hi] in the iv context."""
    return iv.mpf([mpf(lo), mpf(hi)])


def iv_sign(f_iv, lo, hi, dps=120, max_depth=40):
    """Certify a definite sign of f_iv on [lo, hi].

    Returns +1/-1 if every point of a recursive interval cover has that
    definite sign, and None if no definite sign can be certified.
    """
    iv.dps = dps
    return _iv_sign_rec(f_iv, _iv_interval(lo, hi), 0, max_depth)


def _iv_sign_rec(f_iv, u, depth, max_depth):
    try:
        val = f_iv(u)
    except Exception:
        return None
    if val.a > 0:
        return +1
    if val.b < 0:
        return -1
    if depth >= max_depth:
        return None
    a, b = u.a, u.b
    mid = (a + b) / 2
    left = _iv_sign_rec(f_iv, _iv_interval(a, mid), depth + 1, max_depth)
    if left is None:
        return None
    right = _iv_sign_rec(f_iv, _iv_interval(mid, b), depth + 1, max_depth)
    if right is None or right != left:
        return None
    return left


def certify_negative(f_iv, pieces, dps=120, step=None):
    """Certify f_iv < 0 on a list of closed intervals.

    Each piece is split into subintervals of width 'step' (default 1e-4)
    whose signs are certified independently.  Returns True/False and a
    per-piece report.
    """
    iv.dps = dps
    if step is None:
        step = mpf("1e-4")
    ok = True
    report = []
    for lo, hi in pieces:
        lo, hi = mpf(lo), mpf(hi)
        if hi <= lo:
            continue
        n = max(1, int((hi - lo) / step))
        for i in range(n):
            a = lo + (hi - lo) * i / n
            b = lo + (hi - lo) * (i + 1) / n
            s = iv_sign(f_iv, a, b, dps=dps)
            if s != -1:
                ok = False
                report.append((a, b, s))
    return ok, report


def count_zeros(f, lo, hi, N=4000, dps=60):
    """Locate zeros of f on [lo, hi] by a uniform grid + bisection.

    Returns the sorted list of approximate zeros (f(lo), f(hi) itself
    excluded).  Suitable for counting zeros of derivatives on strips.
    """
    mp.dps = dps
    zeros = []
    prev_u, prev_v = None, None
    for i in range(N + 1):
        u = lo + (hi - lo) * i / N
        v = f(u)
        if prev_v is not None and ((v < 0) != (prev_v < 0)):
            lo_b, hi_b = prev_u, u
            for _ in range(200):
                mid = (lo_b + hi_b) / 2
                fm = f(mid)
                if f(lo_b) * fm < 0:
                    hi_b = mid
                else:
                    lo_b = mid
                if hi_b - lo_b < mpf("1e-30"):
                    break
            zeros.append((lo_b + hi_b) / 2)
        prev_u, prev_v = u, v
    return zeros
