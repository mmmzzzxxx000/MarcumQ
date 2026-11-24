# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Validate the manuscript's C(a,nu)*maximum-local-error theorem offline.

Deterministic binary64 computations, exact rational diagnostics, no kernel edits.
The figure tests an absolute-error theorem, not a relative-accuracy guarantee.
"""
import csv
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import random
import sys
import time

HERE = Path(__file__).resolve().parent
from marcum._paths import ROOT as PROJECT
sys.path.insert(0, str(PROJECT / 'src/python'))
from marcum.serr import coefficients

OUT = PROJECT / 'results/section3/recurrence'
NMAX = 256


def require(ok, message):
    if not ok:
        raise ArithmeticError(message)


def fmt(x):
    # Convert only dimensionless ratios or bounded diagnostics to binary64.
    # Tiny absolute errors are stored in scientific notation with mpmath.
    import mpmath as mp
    with mp.workdps(30):
        return mp.nstr(mp.mpf(x.numerator) / x.denominator, 22)


def cases():
    result = [('grid_a%g_v%g' % (a, v), a, v, 'grid')
              for a in (.01, .1, 1., 5., 10., 20., 50., 100.)
              for v in (.1, .25, .5, 1., 2., 10., -.5, -1.5, -2.3)]
    result += [('zero_u1', 2., 2., 'targeted'),
               ('near_zero_u1', math.sqrt(2.), 1., 'targeted'),
               ('near_excluded_order', 1., -1. + 2.**-20, 'targeted')]
    rng = random.Random(20261108)
    for i in range(32):
        a = 10.**rng.uniform(-2, 2)
        v = 10.**rng.uniform(-1, 1) if i % 2 == 0 else -rng.uniform(.05, 4.95)
        result.append(('random_%02d' % i, a, v, 'random'))
    return result


def direct_coefficient(a, v, n):
    # Independent finite Cauchy-product formula, evaluated exactly.
    total = F(0)
    poch = F(1)
    for k in range(n + 1):
        if k:
            poch *= v + k - 1
        total += (F(-1, 2)**(n-k) * (a*a/4)**k /
                  (math.factorial(n-k)*math.factorial(k)*poch))
    return total


def check_case(label, a, nu, group):
    af, vf = F(a), F(nu)
    # A sufficient tail threshold, followed by an exact finite search for m.
    # For n>=safe: denominator>=2n^2 and numerator<=4n+d,
    # so |A_n|+|B_n| <= 2/n+d/(2n^2) <= 1/2.
    d = af*af + 2*abs(vf-2) + 1
    safe = int(math.ceil(max(8., 2*abs(nu-1), math.sqrt(float(2*d))))) + 1
    require(safe >= 8 and F(safe) >= 2*abs(vf-1) and safe*safe >= 2*d,
            label + ': invalid sufficient threshold')
    N = max(NMAX, safe)
    computed = coefficients(a, nu, N)
    require(all(math.isfinite(u) for u in computed), label + ': nonfinite kernel output')
    exact = [F(1), (af*af-2*vf)/(4*vf)]
    H = [F(1), F(1)]
    residual = [F(computed[0])-exact[0], F(computed[1])-exact[1]]
    A, B = [F(0), F(0)], [F(0), F(0)]
    for n in range(2, N+1):
        den = 4*n*(vf+n-1)
        A.append((af*af-2*(vf+2*n-2))/den)
        B.append(1/den)
        exact.append(A[n]*exact[-1]-B[n]*exact[-2])
        residual.append(F(computed[n])-A[n]*F(computed[n-1])+B[n]*F(computed[n-2]))
        H.append(1+abs(A[n])*H[-1]+abs(B[n])*H[-2])
    bad = [n for n in range(2, safe) if abs(A[n])+abs(B[n]) > F(1, 2)]
    m = max(2, max(bad)+1 if bad else 2)
    require(all(abs(A[n])+abs(B[n]) <= F(1, 2) for n in range(m, N+1)),
            label + ': tail contraction failed')
    C = max([F(2)] + H[:m])
    require(all(h <= C for h in H), label + ': H exceeded C')
    for n in (2, 8, 16, 32):
        require(exact[n] == direct_coefficient(af, vf, n), label + ': independent reference failed')
    delta = max(abs(residual[0]), abs(residual[1]))
    rows = []
    first_tiny = None
    max_ratio = F(0)
    max_util = F(0)
    max_normal_util = F(0)
    zeros = 0
    normal_count = 0
    for n in range(2, NMAX+1):
        delta = max(delta, abs(residual[n]))
        error = abs(F(computed[n])-exact[n])
        require(error <= H[n]*delta <= C*delta, label + ': theorem inequality failed')
        tiny = exact[n] != 0 and abs(F(computed[n])) < F(sys.float_info.min)
        if tiny and first_tiny is None:
            first_tiny = n
        ratio = error/delta if delta else F(0)
        util = ratio/H[n]
        max_ratio = max(max_ratio, ratio)
        max_util = max(max_util, util)
        if not tiny:
            normal_count += 1
            max_normal_util = max(max_normal_util, util)
        if error == 0:
            zeros += 1
        rows.append(dict(case=label, group=group, a=a, nu=nu, n=n,
                         computed_u=computed[n], exact_u=fmt(exact[n]),
                         actual_absolute_error=fmt(error), local_residual=fmt(residual[n]),
                         maximum_local_error=fmt(delta), H_n=fmt(H[n]), C_a_nu=fmt(C),
                         absolute_bound_H=fmt(H[n]*delta), absolute_bound_C=fmt(C*delta),
                         error_over_maximum_local=float(ratio), error_over_H_bound=float(util),
                         coefficient_tiny=tiny, after_first_tiny=first_tiny is not None,
                         error_exactly_zero=error == 0, delta_exactly_zero=delta == 0,
                         bound_verified_exactly=True))
    summary = dict(case=label, group=group, a=a, nu=nu, m=m,
                   C_a_nu=fmt(C), H_n_final=fmt(H[NMAX]),
                   checks=len(rows), first_tiny_coefficient=first_tiny,
                   normal_coefficient_checks=normal_count, exact_zero_errors=zeros,
                   max_error_over_local=float(max_ratio), max_error_over_H_bound=float(max_util),
                   max_normal_error_over_H_bound=float(max_normal_util),
                   violations=0, direct_formula_checks=4)
    return rows, summary


def write_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    kernel = PROJECT/'src/python/marcum/serr.py'
    before = hashlib.sha256(kernel.read_bytes()).hexdigest()
    started = time.perf_counter()
    all_rows, summaries = [], []
    for label, a, v, group in cases():
        rows, summary = check_case(label, a, v, group)
        all_rows.extend(rows)
        summaries.append(summary)
        print('%s: %d exact checks, max error/H-bound=%.5g, C=%s' %
              (label, len(rows), summary['max_error_over_H_bound'], summary['C_a_nu']), flush=True)
    require(before == hashlib.sha256(kernel.read_bytes()).hexdigest(), 'Kernel changed')
    write_csv(OUT/'steps.csv', all_rows)
    write_csv(OUT/'cases.csv', summaries)
    summary = dict(cases=len(summaries), checks=len(all_rows), Nmax=NMAX,
                   violations=0, reference='Exact rational recurrence plus independent finite Cauchy product',
                   input_interpretation='Exact values of binary64 input parameters',
                   residual_arithmetic='Exact rational, including initialization errors',
                   bound_arithmetic='Exact rational; decimal strings are display only',
                   includes_underflow=True, underflow_marked_separately=True,
                   checks_before_first_tiny=sum(not r['after_first_tiny'] for r in all_rows),
                   checks_after_first_tiny=sum(r['after_first_tiny'] for r in all_rows),
                   maximum_error_over_H_bound=max(s['max_error_over_H_bound'] for s in summaries),
                   random_seed=20261108, kernel_sha256=before,
                   elapsed_seconds=time.perf_counter()-started)
    (OUT/'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
