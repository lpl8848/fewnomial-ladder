# fewnomial-ladder

Companion code for the paper

> **The fewnomial ladder: $2t-1$ positive solutions for trinomial systems,
> and the sharpness problem for $\mathcal{S}_+(3,t)$** — Peilin Luo.

The code reproduces every numerical and exact certificate used in the paper:

| script | what it certifies | paper reference |
|---|---|---|
| `scripts/certify_base.py` | the five roots of $F_1$ in $(0,1)$, each in a certified bracket of width $\approx 2\cdot10^{-95}$ with opposite interval signs at 120 digits; the certified values $F_1(0.01)$, $F_1(0.55)$, $F_1(0.5761)$; the four cluster extrema and the margin $m = 2.6\cdot10^{-4}$; the six marks $e_0,\dots,e_5$ | Lemma *base* |
| `scripts/certify_instances.py` | $F_4, F_5, F_6$ have exactly $7, 9, 11$ roots in $(0,1)$: certified brackets (existence), pointwise negativity on $K = [0.01,0.55]\setminus\bigcup_j I_j$ and on $[0,0.01]$ (exhaustiveness), mark sign preservation, cluster and bump-strip uniqueness | Theorem *verified* |
| `scripts/certify_m3.py` | **exact** certificate $M_3\ge 5$: the resultant $R(y)=P(y^6)^3-y^3(1-y^6)^{49}Q(y^6)^3$ of degree 702 (its square is the resultant of the system in $z$), Sturm count of exactly five roots in $(0,1)$, and membership of the five $\psi$-roots | Theorem $M_3\ge5$ |
| `scripts/certify_table.py` | Table 1 (bump peaks, amplitudes, values at $0.5761$) and the mark perturbations of the pointwise certification | Table 1 |
| `scripts/verify_reduction.py` | symbolic check that the Müller–Regensburger system reduces exactly to $F_1(u)=0$ with $u=bx$, $v=y^5/x$ | Lemma *reduction* |

The two general lemmas of the paper (single-bump localization and
multi-bump localization) are pure existence statements (take the exponents
large enough); they are proved in the paper and need no code.

## Requirements

- Python 3.9+ (`py -3.9` on Windows)
- `mpmath >= 1.3`
- `sympy >= 1.12`

Install with `pip install -r requirements.txt`.

## Running

From the repository root:

```
python scripts/certify_base.py
python scripts/certify_instances.py
python scripts/certify_m3.py
python scripts/certify_table.py
python scripts/verify_reduction.py
```

All scripts end with an explicit `CHECK PASSED` statement (and `assert`
every certificate), so a clean run is itself the certificate.

## Precision and rigor — read this

- Root brackets are refined by bisection to width $\approx 2\cdot10^{-95}$
  and their signs are re-verified with mpmath interval arithmetic at 120
  digits. mpmath's interval context is **midpoint–radius** arithmetic and
  does not perform outward rounding; on the thin brackets used here every
  sign is separated from zero by a margin far above the evaluation error.
- The count of Theorem $M_3\ge5$ does **not** rely on floating point: it is
  certified exactly by a sympy resultant + Sturm computation over
  $\mathbb{Q}$ (`scripts/certify_m3.py`).
- For a formally outward-rounded version of the interval certificates, swap
  the interval backend for `arb`/`MPFI`; the statements and outputs are
  unchanged.

## Layout

```
fewnomial_ladder/
    functions.py    # F1, bumps, F4/F5/F6, psi (+ interval-valued variants),
                    # certified marks, margins and localization intervals
    certify.py      # logit-grid discovery, bisection, interval sign covers
scripts/            # the five runnable certificates listed above
```

## License

MIT — see `LICENSE`.
