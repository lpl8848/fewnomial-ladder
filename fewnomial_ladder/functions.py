# -*- coding: utf-8 -*-
"""Core functions of ``The fewnomial ladder`` paper, evaluated with mpmath.

All constants match ``fewnomial_ladder.tex``:

* ``b = 1.392 = 174/125``;
* ``F1`` (eq. (F1)), the one-variable form of the Mueller--Regensburger system;
* the monomial bumps ``amp * u**p * (1-u)**q`` of eqs. (F4)--(F6);
* ``psi`` (eq. (m3ex)), the five-root example for the M_3 family;
* the certified marks ``e_0..e_5`` and the cluster margin ``m`` of Lemma ``base``;
* the localization intervals used in the certification of F_4, F_5, F_6.

Point-valued functions take/return mpmath scalars (``mp.dps`` context);
``*_iv`` variants take/return mpmath intervals (``iv.dps`` context).
"""

from mpmath import mpf


def _b():
    """b = 1.392 = 174/125, evaluated at the current mp context precision."""
    return mpf(174) / 125


# ----------------------------------------------------------------------
# Point-valued functions (mp context)
# ----------------------------------------------------------------------
def F1(u):
    """One-variable form of the Mueller--Regensburger system (eq. (F1))."""
    b = _b()
    return (b ** (-mpf(24) / 5) * u ** (mpf(24) / 5) * (1 - u) ** (-mpf(1) / 5)
            + b ** (mpf(4) / 5) * u ** (mpf(1) / 5) * (1 - u) ** (mpf(1) / 5) - 1)


def F1p(u):
    """Derivative of F1, written explicitly (used with mp.diff elsewhere too)."""
    from mpmath import diff
    return diff(F1, u)


def bump(amp, p, q, u):
    """Monomial bump amp * u^p * (1-u)^q, p, q nonnegative integers."""
    return amp * u ** p * (1 - u) ** q


def bumpp(amp, p, q, u):
    """Derivative of the bump: B'(u) = B(u) * (p/u - q/(1-u))."""
    return bump(amp, p, q, u) * (mpf(p) / u - mpf(q) / (1 - u))


# Bump parameters (amp, p, q, u_*) of eqs. (F4)--(F6); u_* = p/(p+q).
BUMP_DEFS = {
    4: [(10 ** 15, 15, 45)],
    5: [(1400, 2, 44), (6 * 10 ** 7, 8, 32)],
    6: [(1400, 2, 44), (6 * 10 ** 7, 8, 32), (4 * 10 ** 33, 50, 69)],
}


def Ft(t, u):
    """F_t(u) = F1(u) + sum of the bumps of the t-th instance (t = 4, 5, 6)."""
    total = F1(u)
    for amp, p, q in BUMP_DEFS[t]:
        total += bump(amp, p, q, u)
    return total


def psi(x):
    """The M_3 example (eq. (m3ex)): x^{-1/6}(1-x)^{-49/3} P(x) - Q(x)."""
    P = x ** 3 - mpf(617) / 500 * x ** 2 + mpf(441) / 1000 * x - mpf(43) / 1000
    Q = x ** 3 - mpf(421) / 250 * x ** 2 + mpf(413) / 500 * x - mpf(107) / 1000
    return x ** (-mpf(1) / 6) * (1 - x) ** (-mpf(49) / 3) * P - Q


# ----------------------------------------------------------------------
# Interval-valued functions (iv context)
# ----------------------------------------------------------------------
def _iv_b():
    from mpmath import iv
    return iv.mpf(174) / 125


def F1_iv(u):
    from mpmath import iv
    b = _iv_b()
    return (b ** (-iv.mpf(24) / 5) * u ** (iv.mpf(24) / 5) * (1 - u) ** (-iv.mpf(1) / 5)
            + b ** (iv.mpf(4) / 5) * u ** (iv.mpf(1) / 5) * (1 - u) ** (iv.mpf(1) / 5) - 1)


def F1p_iv(u):
    """d/du F1: both terms differentiated explicitly."""
    from mpmath import iv
    b = _iv_b()
    t1 = b ** (-iv.mpf(24) / 5) * u ** (iv.mpf(19) / 5) * (1 - u) ** (-iv.mpf(6) / 5) \
        * (iv.mpf(24) / 5 * (1 - u) + iv.mpf(1) / 5 * u)
    t2 = b ** (iv.mpf(4) / 5) * u ** (-iv.mpf(4) / 5) * (1 - u) ** (-iv.mpf(4) / 5) \
        * (iv.mpf(1) / 5 * (1 - u) - iv.mpf(1) / 5 * u)
    return t1 + t2


def bump_iv(amp, p, q, u):
    return iv_of(amp) * u ** p * (1 - u) ** q


def bumpp_iv(amp, p, q, u):
    from mpmath import iv
    return bump_iv(amp, p, q, u) * (iv.mpf(p) / u - iv.mpf(q) / (1 - u))


def iv_of(x):
    """Convert an exact rational/integer amplitude into the iv context."""
    from mpmath import iv
    from fractions import Fraction
    if isinstance(x, int):
        return iv.mpf(x)
    f = Fraction(str(x)) if not isinstance(x, Fraction) else x
    return iv.mpf(f.numerator) / f.denominator


def Ft_iv(t, u):
    total = F1_iv(u)
    for amp, p, q in BUMP_DEFS[t]:
        total += bump_iv(amp, p, q, u)
    return total


def Ftp_iv(t, u):
    total = F1p_iv(u)
    for amp, p, q in BUMP_DEFS[t]:
        total += bumpp_iv(amp, p, q, u)
    return total


# ----------------------------------------------------------------------
# Certified constants of the paper
# ----------------------------------------------------------------------
# Marks e_0..e_5 of Lemma "base" (cluster, with alternating certified signs).
MARKS = [mpf("0.55"), mpf("0.6435"), mpf("0.8041"), mpf("0.9035"),
         mpf("0.9712"), mpf("0.99")]

# Cluster margin of Lemma "base" (c): every inter-root extremum has |F1| >= m.
MARGIN_M = mpf("2.6e-4")

# Localization intervals I_j (open) used in the certification of the
# instances, pairwise disjoint and inside (0, 0.55); each contains the two
# roots created by the corresponding bump.  For F_5/F_6 bump (i) the left
# endpoint sits below 0.01, so the region [0, 0.01] is certified separately.
INTERVALS = {
    4: [(mpf("0.03"), mpf("0.47"))],
    5: [(mpf("0.0015"), mpf("0.0855")), (mpf("0.1"), mpf("0.3"))],
    6: [(mpf("0.0015"), mpf("0.0855")), (mpf("0.1"), mpf("0.3")),
        (mpf("0.30017"), mpf("0.54017"))],
}

# K-region of the certification: K = [0.01, 0.55] \ union_j I_j, written as a
# list of closed subintervals (the two endpoints of I_j themselves carry
# certified signs and are included in the cover).
K_LOW, K_HIGH = mpf("0.01"), mpf("0.55")


def k_regions(t):
    """Closed subintervals covering [0.01, 0.55] minus the open I_j's."""
    pieces = [(K_LOW, K_HIGH)]
    for lo, hi in INTERVALS[t]:
        new = []
        for a, b in pieces:
            if hi <= a or lo >= b:
                new.append((a, b))
            else:
                if a < lo:
                    new.append((a, lo))
                if hi < b:
                    new.append((hi, b))
        pieces = new
    return pieces
