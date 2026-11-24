# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Direct cancellation diagnostics for unchanged Algorithm 1 at fixed N."""
import hashlib
import json
from pathlib import Path
import sys

import mpmath as mp

sys.dont_write_bytecode = True
from marcum._paths import ROOT
OUT = ROOT / "results/section4/cancellation"
OUT.mkdir(parents=True, exist_ok=True)
from marcum._paths import ROOT as PROJECT
sys.path.insert(0, str(PROJECT/'src/python'))
from marcum.serr import evaluate as binary64

N = 128


def exact_quantities(a0, b0, v0, dps):
    with mp.workdps(dps):
        a, b, v = map(mp.mpf, (a0, b0, v0))
        u = [mp.mpf(1), (a*a-2*v)/(4*v)]
        for n in range(2, N+1):
            den = 4*n*(v+n-1)
            u.append((a*a-2*(v+2*n-2))/den*u[-1]-u[-2]/den)
        prefactor = mp.exp(-a*a/2)/(2**(v-1)*mp.gamma(v))
        terms = [prefactor*u[n]*b**(2*(v+n))/(2*(v+n)) for n in range(N+1)]
        total = mp.fsum(terms)
        mass = mp.fsum(abs(t) for t in terms)
        q = 1-total
        m = 2*mp.exp(-a*a/4+b*b/2)*(b*b/2)**v/mp.gamma(v+1)
        upper = ((1+m)/(1-m) if m < mp.mpf('.998') else
                 2*(1+m)*mp.gamma(v)/mp.gammainc(v, b*b/2, mp.inf))
        tail = (2*mp.exp(-a*a/4)*(b*b/2)**v/mp.gamma(v+1)
                *(b*b/2)**(N+1)/mp.factorial(N+1)/(1-(b*b/2)/(N+2)))
        return dict(a=a, b=b, nu=v, N=N, M=m, B=upper,
                    sum_abs=mass, sum=total, Q_N=q, kappa=(1+mass)/abs(q),
                    positive_sum=mp.fsum(t for t in terms if t > 0),
                    negative_magnitude=mp.fsum(-t for t in terms if t < 0),
                    max_abs_term=max(abs(t) for t in terms), tail_bound=tail)


def main():
    source = PROJECT/'src/python'/'marcum'/'serr.py'
    before = hashlib.sha256(source.read_bytes()).hexdigest()
    result = []
    for case in ((8, 5, 1), (3, 2.556, 1), (1, 7, 1)):
        low = exact_quantities(*case, dps=100)
        high = exact_quantities(*case, dps=180)
        with mp.workdps(180):
            mismatch = max(abs(low[key]/high[key]-1) for key in ('sum_abs', 'Q_N', 'kappa'))
            assert mismatch < mp.mpf('1e-75')
            assert high['B'] >= high['kappa']
            assert high['tail_bound']/abs(high['Q_N']) < mp.mpf('1e-25')
            actual = binary64(*map(float, case), N)
            high['Q_double'] = mp.mpf(actual)
            high['float_relative_error'] = abs(mp.mpf(actual)-high['Q_N'])/abs(high['Q_N'])
            high['precision_difference'] = mismatch
            a, b, v = map(mp.mpf, case)
            def density(t):
                return a**(1-v)*t**v*mp.exp(-(t*t+a*a)/2)*mp.besseli(v-1, a*t)
            with mp.workdps(80):
                integral = mp.quad(density, [b, b+2, max(a, b)+12, mp.inf])
            high['Bessel_reference'] = integral
            high['truncation_relative_error'] = abs(high['Q_N']-integral)/abs(integral)
            assert high['truncation_relative_error'] < mp.mpf('1e-25')
            result.append({key: (mp.nstr(value, 55) if isinstance(value, mp.mpf) else value)
                           for key, value in high.items()})
            print(json.dumps(result[-1], indent=2), flush=True)
    assert hashlib.sha256(source.read_bytes()).hexdigest() == before
    (OUT/'results.json').write_text(json.dumps(dict(
        method='Unchanged Python Algorithm 1 in binary64; fixed N=128 in both cases.',
        source_sha256=before, high_precision_digits=[100, 180],
        independent_Bessel_integral_digits=80, cases=result), indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
