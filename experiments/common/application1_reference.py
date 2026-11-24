# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Neutral direct lower Bessel integral at every Application 1 threshold.

No SE, Poisson mixture, incomplete-Gamma sum, or SE-RR used as oracle.
All manuscript cases have 2*nu-1 integral and >=1, so the density is analytic
at t=0. Independent quadrature orders/precisions and tanh-sinh spot checks.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import mpmath as mp

PROTOCOL = 'application1_direct_lower_bessel_50_80dps_gl8_12_v1'


@lru_cache(maxsize=4)
def rule(n, dps):
    with mp.workdps(dps):
        xx, ww = mp.gauss_quadrature(n, 'legendre')
        return tuple(xx), tuple(ww)


def distribution(case, bs):
    versions = []
    for dps, nodes in ((50, 8), (80, 12)):
        with mp.workdps(dps):
            a, nu = mp.mpf(case['a']), mp.mpf(case['nu'])
            if nu <= 0 or int(2*nu-1) != 2*nu-1 or nu < 1:
                raise ValueError('This reference quadrature is validated for the six manuscript cases only')
            xx, ww = rule(nodes, dps)
            def density(t):
                return a**(1-nu)*t**nu*mp.exp(-(t*t+a*a)/2)*mp.besseli(nu-1, a*t)
            total, left = mp.mpf(0), mp.mpf(0)
            vals = []
            for b in bs:
                right = mp.mpf(float(b))
                if right > left:
                    mid, half = (left+right)/2, (right-left)/2
                    total += half*mp.fsum(w*density(mid+half*x) for x, w in zip(xx, ww))
                vals.append(+total)
                left = right
            versions.append(vals)
    with mp.workdps(100):
        differences = [abs(v-w)/abs(w) for v, w in zip(*versions) if w]
        maximum = max(differences) if differences else mp.mpf(0)
        if maximum > mp.mpf('1e-25'):
            raise ArithmeticError('CDF integral configurations disagree: '+str(maximum))
        a, nu = mp.mpf(case['a']), mp.mpf(case['nu'])
        def density(t):
            return a**(1-nu)*t**nu*mp.exp(-(t*t+a*a)/2)*mp.besseli(nu-1, a*t)
        checks = []
        for i in sorted(set((1, len(bs)//2, len(bs)-1))):
            b = mp.mpf(float(bs[i]))
            independent = mp.quad(density, [0, b/2, b], method='tanh-sinh')
            diff = abs(independent-versions[1][i])/abs(independent) if independent else mp.mpf(0)
            if diff > mp.mpf('1e-35'):
                raise ArithmeticError('Independent CDF tanh-sinh check disagrees')
            checks.append(dict(index=i, relative_difference=mp.nstr(diff, 30)))
        return dict(protocol=PROTOCOL, case=case, b=[float(b) for b in bs],
                    p=[mp.nstr(v, 78) for v in versions[1]],
                    low_p=[mp.nstr(v, 48) for v in versions[0]],
                    configurations=[dict(dps=50, gauss_nodes_per_interval=8), dict(dps=80, gauss_nodes_per_interval=12)],
                    max_relative_configuration_difference=mp.nstr(maximum, 30),
                    independent_tanh_sinh_checks=checks,
                    uses_se_or_serr=False, interval_certificate=False)


def run(folder, cases, bs, workers=3):
    folder = Path(folder)/'references'
    folder.mkdir(parents=True, exist_ok=True)
    fingerprint = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    found, pending = {}, []
    for case in cases:
        path = folder/(case['id']+'.json')
        if path.is_file():
            old = json.loads(path.read_text(encoding='utf-8'))
            if old.get('source_sha256') == fingerprint and old['case'] == case and old['b'] == list(bs):
                found[case['id']] = old
                continue
        pending.append(case)
    if pending:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            jobs = [pool.submit(distribution, case, list(bs)) for case in pending]
            for job in as_completed(jobs):
                record = job.result()
                record['source_sha256'] = fingerprint
                key = record['case']['id']
                (folder/(key+'.json')).write_text(json.dumps(record, indent=2), encoding='utf-8')
                found[key] = record
                print('Neutral direct CDF integral:', key, 'complete', flush=True)
    return found


if __name__ == '__main__':
    import argparse
    import numpy as np
    root = Path(__file__).resolve().parents[1]
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', default=str(root/'results/application1/rebuilt_20261030'))
    p.add_argument('--workers', type=int, default=3)
    args = p.parse_args()
    cfg = json.loads((root/'configs/application1_rebuilt.json').read_text(encoding='utf-8'))
    run(args.output, cfg['cases'], np.linspace(*cfg['b_range'], cfg['points']), args.workers)
