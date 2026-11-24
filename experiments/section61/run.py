# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Fortran CDF timings with GammaCHI and hybrid settings selected to 1e-10."""

import argparse
import json
from pathlib import Path

from helpers import (APP, OUT, ROOT, ct, exact_error, hybrid_mask,
                     math, mp, np, power_status, random, sha, shutil, statistics,
                     subprocess, sys, write_csv, write_json)
from application1_reference import run as reference_grid

sys.path.insert(0, str(ROOT / 'experiments/common'))
from build_rebuilt import UPSTREAM, fortran_text, replace_once

MODES = ('scalar', 'batch')


def gamma_source(original):
    text = original.read_text(encoding='utf-8')
    upstream = UPSTREAM / '939/Fortran90/Src/marcumq.f90'
    _, changes = fortran_text(upstream.read_text(encoding='utf-8'))
    edits = []
    for item in changes[1:5]:
        text = replace_once(text, item['original'], item['replacement'], edits)
    cache = (ROOT / 'src/fortran/application1_cache.f90').read_text(encoding='utf-8')
    declarations, _, procedures = cache.partition('\ncontains\n')
    declarations = declarations[declarations.index(' logical ::'):]
    procedures = procedures[:procedures.rindex('\nend module')]
    text = replace_once(text, '             gamstar, quotgamm\n  CONTAINS',
                        '             gamstar, quotgamm\n'
                        '  public :: qser, pser, bench_setting, app1_cached, app1_cache_begin\n'
                        '  integer :: bench_setting=-1\n'+declarations+
                        '\n  CONTAINS\n'+procedures+'\n', edits)
    for name in ('pser', 'qser'):
        start = text.index('    SUBROUTINE '+name+'(')
        stop = text.index('    END SUBROUTINE '+name, start)
        body = text[start:stop]
        for old, new in [('loggam(mu+1.0_r8)', 'cdf_lgmu(mu,0)'),
                         ('loggam(mu+n+1.0_r8)', 'cdf_lgmu(mu,n)'),
                         ('loggam(n+1.0_r8)', 'cdf_lgn(n)'),
                         ('log(x)', 'cdf_logx(x)'), ('exp(-x)', 'cdf_expx(x)')]:
            if old in body:
                edits.append(dict(scope=name, original=old, replacement=new, occurrences=body.count(old)))
                body = body.replace(old, new)
        text = text[:start]+body+text[stop:]
    return text, edits


def build(method, output):
    dest = output / method
    dest.mkdir(parents=True, exist_ok=True)
    compiler = shutil.which('gfortran')
    options = ['-O3', '-frecursive', '-fno-fast-math', '-ffp-contract=off', '-fPIC',
               '-ffree-line-length-none', '-shared', '-static-libgfortran', '-static-libgcc']
    if method == 'GammaCHI':
        original = ROOT / 'third_party/gammachi/gammaCHI.f90'
        constants = ROOT / 'third_party/gammachi/cons.f90'
        generated = dest / 'gammachi_generated.f90'
        text, edits = gamma_source(original)
        generated.write_text(text, encoding='utf-8')
        sources = [constants, generated, OUT / 'gammachi_fixed.f90']
        options += ['-I', str(original.parent)]
    else:
        sources = [APP / (name+'.o') for name in
                   ('cons', 'incgam', 'application1_cache', 'marcumq_application1', 'application1')]
        sources.append(OUT / 'hybrid_fixed.f90')
        options += ['-I', str(APP)]
        edits = []
    path = dest / 'fixed.dll'
    command = [compiler]+options+[str(source) for source in sources]+['-o', str(path)]
    subprocess.run(command, cwd=dest, check=True)
    lib = ct.CDLL(str(path))
    double = np.ctypeslib.ndpointer(np.float64, ndim=1, flags='C_CONTIGUOUS')
    int32 = np.ctypeslib.ndpointer(np.int32, ndim=1, flags='C_CONTIGUOUS')
    int64 = np.ctypeslib.ndpointer(np.int64, ndim=1, flags='C_CONTIGUOUS')
    sig = [ct.c_int, ct.c_double, ct.c_double, ct.c_int64, double, int32, int64, double, int32]
    lib.fixed_values.argtypes, lib.fixed_values.restype = sig, None
    lib.fixed_time.argtypes, lib.fixed_time.restype = sig+[ct.c_int64], ct.c_double
    return lib, dict(command=command, edits=edits, sources={str(s): sha(s) for s in sources},
                     library_sha256=sha(path))


def prepare(lib, mode, case, b, branches, orders):
    values = np.empty(len(b))
    flags = np.empty(len(b), dtype=np.int32)
    args = (int(mode == 'batch'), case['a'], case['nu'], len(b), b, branches, orders, values, flags)
    def evaluate():
        lib.fixed_values(*args)
        return values.copy(), flags.copy()
    return evaluate, lambda n: lib.fixed_time(*args, n)


def calibrate(lib, mode, case, b, branches, reference, target, first=None):
    orders = np.full(len(b), -1, dtype=np.int64)
    if first is None:
        first = np.zeros(len(b), dtype=np.int64)
    for order in range(1001):
        evaluate, _ = prepare(lib, mode, case, b, branches, np.maximum(first, order))
        values, flags = evaluate()
        for i in np.flatnonzero((orders < 0) & (flags == 0) & (first <= order)):
            if exact_error(values[i], reference[i]) <= target:
                orders[i] = order
        if np.all(orders >= 0):
            return orders
    raise RuntimeError('Calibration failed for '+case['id'])


def main(grid_max=5):
    output = ROOT / 'results/section61'
    output.mkdir(parents=True, exist_ok=True)
    cfg = __import__('json').loads((ROOT/'configs/section61/cdf.json').read_text(encoding='utf-8'))
    b = np.linspace(0, grid_max, cfg['points'])
    stops = (('0_to_'+format(grid_max, 'g'), len(b)),)
    reference_grid(output, cfg['cases'], b, workers=3)
    ref_dir=output/'references'
    reference_hashes={case['id']:sha(ref_dir/(case['id']+'.json')) for case in cfg['cases']}
    methods = ('GammaCHI', 'Hybrid')
    dlls, manifests = {}, {}
    for method in methods:
        dlls[method], manifests[method] = build(method, output)
    jobs = []
    accuracy = []
    with mp.workdps(100):
        target = mp.mpf(cfg['target_rtol'])
        for case in cfg['cases']:
            path = ref_dir / (case['id']+'.json')
            assert sha(path) == reference_hashes[case['id']]
            reference_record = json.loads(path.read_text(encoding='utf-8'))
            assert reference_record['b'] == list(b)
            ref = [mp.mpf(p) for p in reference_record['p']]
            x = case['a']**2/2
            assert x < 30 and case['a'] <= 7
            tails = (b*b/2 > x+case['nu']).astype(np.int32)
            route = hybrid_mask(case['a'], case['nu'], b).astype(np.int32)
            for mode in MODES:
                gamma_orders = calibrate(dlls['GammaCHI'], mode, case, b, tails, ref, target)
                hybrid = calibrate(dlls['Hybrid'], mode, case, b, route, ref, target,
                                   first=route.astype(np.int64))
                for method, branches, orders in (('GammaCHI', tails, gamma_orders), ('Hybrid', route, hybrid)):
                    for interval, count in stops:
                        evaluate, timer = prepare(dlls[method], mode, case, b[:count], branches[:count], orders[:count])
                        values, flags = evaluate()
                        errors = [exact_error(v, truth) for v, truth in zip(values, ref)]
                        assert np.all(flags == 0) and max(errors) <= target, (case, mode, method, interval)
                        row = dict(case=case['id'], M=case['nu'], a=case['a'], interval=interval,
                                   points=count, se_rr_points=int(route[:count].sum()), method=method,
                                   mode=mode, max_relative_error=float(max(errors)), samples=[])
                        jobs.append((row, timer))
                    accuracy.append(dict(case=case['id'], method=method, mode=mode, orders=orders.tolist()))
            print('Calibrated '+case['id'], flush=True)
    power_before = power_status()
    for row, timer in jobs:
        row['repeats'] = max(4, math.ceil(.016/timer(16)))
    rng = random.Random(20261111)
    for _ in range(9):
        rng.shuffle(jobs)
        for row, timer in jobs:
            row['samples'].append(1e6*timer(row['repeats'])/row['points'])
    power_after = power_status()
    for row, _ in jobs:
        row['microseconds_per_value'] = statistics.median(row['samples'])
    timings = {(r['case'], r['interval'], r['mode'], r['method']): r for r, _ in jobs}
    pairs = []
    for case in cfg['cases']:
        for interval, _ in stops:
            for mode in MODES:
                base, hybrid = [timings[case['id'], interval, mode, method] for method in methods]
                t1, t2 = (r['microseconds_per_value'] for r in (base, hybrid))
                pairs.append(dict(case=case['id'], M=case['nu'], a=case['a'], interval=interval, mode=mode,
                                  points=base['points'], se_rr_points=base['se_rr_points'],
                                  gammachi_us=t1, hybrid_us=t2, speedup=t1/t2,
                                  gammachi_max_error=base['max_relative_error'],
                                  hybrid_max_error=hybrid['max_relative_error']))
    write_csv(output / 'pairs.csv', pairs)
    write_json(output / 'summary.json', dict(pairs=pairs, trials=[r for r, _ in jobs], orders=accuracy,
               builds=manifests, power_before=power_before, power_after=power_after,
               target=cfg['target_rtol'], reference_sha256=reference_hashes,
               b_range=[float(b[0]), float(b[-1])], grid_points=len(b),
               protocol='Nine interleaved Fortran medians. Parameters, branches and orders set outside timers.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--grid-max', type=float, default=5)
    main(parser.parse_args().grid_max)
