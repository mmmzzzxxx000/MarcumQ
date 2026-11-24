# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Time the methods selected by the two hybrid strategies on the Section 5.1 grid."""
from collections import Counter
import csv
import ctypes
import functools
import hashlib
import json
import math
from pathlib import Path
import random
import statistics
import subprocess
import sys

import mpmath as mp
import numpy as np


sys.dont_write_bytecode = True
OUT = Path(__file__).resolve().parent
from marcum._paths import ROOT as PROJECT
SOURCE = OUT
OUT = PROJECT / 'results/section52/grid'
OUT.mkdir(parents=True, exist_ok=True)
BENCHMARK = PROJECT / 'data/section51'
sys.path[:0] = [str(PROJECT / part) for part in ('src', 'experiments', 'validation')]
sys.path.insert(0, str(OUT.parent / 'B_bound_triptych_20261109'))
from common import power_status
from hybrid_runtime import load as load_fortran
from marcum.rebuilt import kernel as original_kernel
from fp_environment import load as load_fp_guard
from bound import upper
from marcum.rr import kernel as python_rr_kernel


LANGUAGES = ('python', 'cpp', 'fortran')
ORDERS = (1, 10, 20)
STRATEGIES = ('algorithm_939', 'proposed')
METHODS = ('se', 'af', 'ni', 'se_rr', 'rr')
GRID = np.arange(1, 41) / 2


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_cpp_rr():
    native = PROJECT/'build/rebuilt/cpp/components_generated.cpp'
    source = OUT/'components_rr.cpp'
    source.write_text(native.read_text(encoding='utf-8')+'\n#include "rr.inc"\n', encoding='utf-8')
    library = OUT/'rr_port.dll'
    command = ['g++', '-O3', '-fno-fast-math',
               '-ffp-contract=off', '-std=c++17', '-shared', str(source),
               '-I', str(PROJECT/'src/cpp'), '-I', str(SOURCE), '-I', str(PROJECT/'src/cpp'),
               str(PROJECT/'src/cpp/serr.cpp'), '-static-libgcc',
               '-static-libstdc++', '-o', str(library)]
    subprocess.run(command, cwd=str(OUT), check=True)
    dll = ctypes.CDLL(str(library))
    signature = [ctypes.c_double]*3+[ctypes.c_int64]
    dll.rr_port_q.argtypes = signature
    dll.rr_port_q.restype = ctypes.c_double
    dll.rr_port_time.argtypes = signature+[ctypes.c_int64]
    dll.rr_port_time.restype = ctypes.c_double
    return dll, command


def plans(fortran, a, b, nu):
    x, y = a*a/2, b*b/2
    w = np.sqrt(4*x+2*nu)
    af = (a*b > 30) & (nu*nu < 2*a*b)
    baseline = np.select([x<30, af, (y>=x+nu-w)&(y<=x+nu+w)], [0, 1, 4], default=2)
    bound, _ = upper(a, b, nu)
    proposed = np.select([(b<=7)&(bound<=1000), (b<=7)|(a<=7), af,
                          (y>x+nu-w)&(y<x+nu+w)], [3, 0, 1, 4], default=2)
    diagonal = (a==b)
    for plan in (baseline, proposed):
        plan[(plan==1)&diagonal] = 0
    return baseline, proposed


def read_calibration():
    audit = json.loads((BENCHMARK/'audit.json').read_text(encoding='utf-8'))
    assert audit['complete'] and not audit['issues']
    path = BENCHMARK/'calibration.jsonl'
    assert sha(path)==audit['artifact_sha256'][path.name]
    with path.open(encoding='utf-8') as source:
        entries = [json.loads(line) for line in source]
    return {(row['language'], row['method'], row['a'], row['b'], row['nu']): row
            for row in entries}


def rr_limit(a, b, nu):
    x, y = a*a/2, b*b/2
    qser_fallback = y>x+nu and y-x+(1-math.sqrt(2*(x+y)+1))<5 and x<200
    return 1000 if qser_fallback else 10


def main():
    fp_guard, _ = load_fp_guard()
    fortran, _ = load_fortran('rr')
    cpp, cpp_command = build_cpp_rr()
    reference_calibration = read_calibration()
    fp_guard.million_set_binary64()
    a, b = np.meshgrid(GRID, GRID)
    rows = []
    for nu in ORDERS:
        for strategy, plan in zip(STRATEGIES, plans(fortran, a, b, nu)):
            for row, col in np.ndindex(plan.shape):
                method = METHODS[plan[row,col]]
                for language in LANGUAGES:
                    rows.append(dict(language=language, strategy=strategy, nu=nu,
                                     a=float(a[row,col]), b=float(b[row,col]), method=method))

    @functools.lru_cache(maxsize=None)
    def calibrate(language, method, a, b, nu):
        reference = mp.mpf(reference_calibration[language,'se',a,b,nu]['reference'])
        if method != 'rr':
            evaluate, timer = original_kernel(method,a,b,float(nu),language)
            setting = reference_calibration[language,method,a,b,nu]['setting']
            if setting is None:
                return setting, None, timer
        else:
            if language == 'python':
                evaluate, timer = python_rr_kernel(a,b,float(nu))
            elif language == 'cpp':
                evaluate = lambda k: cpp.rr_port_q(a,b,float(nu),k)
                timer = lambda k,r: cpp.rr_port_time(a,b,float(nu),k,r)
            else:
                evaluate = lambda k: fortran.rebuilt_q(4,a,b,float(nu),k)
                timer = lambda k,r: fortran.rebuilt_time(4,a,b,float(nu),k,r)
            setting = None
            for trial in range(rr_limit(a,b,nu)+1):
                fp_guard.million_set_binary64()
                value = evaluate(trial)
                if math.isfinite(value) and abs(mp.mpf(value)-reference)<=mp.mpf('1e-10')*reference:
                    setting = trial
                    if language != 'fortran':
                        native = fortran.rebuilt_q(4,a,b,float(nu),trial)
                        assert math.isfinite(native) and abs(value-native)<1e-12, (a,b,nu,trial,value,native)
                    break
            if setting is None:
                return setting, None, timer
        fp_guard.million_set_binary64()
        error = abs(mp.mpf(evaluate(setting))-reference)/reference
        assert mp.isfinite(error) and error<=mp.mpf('1e-10')
        return setting, float(error), timer

    def workload(row):
        return tuple(row[k] for k in ('language','method','a','b','nu'))

    jobs = list(dict.fromkeys(workload(row) for row in rows))
    random.Random(20261109).shuffle(jobs)
    with mp.workdps(100):
        calibrated = {job: calibrate(*job) for job in jobs}
    failures = [job for job, (setting, _, _) in calibrated.items() if setting is None]
    print('Calibrated %d unique workloads for %d cells; %d failures.' % (len(jobs),len(rows),len(failures)), flush=True)
    if failures:
        print('Failures:', failures[:30], flush=True)
    power_before = power_status()
    results = {}
    for index, job in enumerate(jobs):
        setting, error, timer = calibrated[job]
        if setting is None:
            results[job] = dict(setting=None, relative_error=None, microseconds=None, samples=[], repeats=0)
            continue
        fp_guard.million_set_binary64()
        pilot = timer(setting,10)
        repeats = max(10,min(30000,math.ceil(.0006/max(pilot,1e-9))))
        samples = []
        for _ in range(5):
            fp_guard.million_set_binary64()
            elapsed = timer(setting,repeats)*1e6
            assert math.isfinite(elapsed) and elapsed>0
            samples.append(elapsed)
        results[job] = dict(setting=setting, relative_error=error,
                            microseconds=statistics.median(samples),samples=samples,repeats=repeats)
        if (index+1)%2000==0:
            print('Timed %d/%d' % (index+1,len(jobs)), flush=True)
    power_after = power_status()
    for row in rows:
        row.update(results[workload(row)])
    assert len(rows)==28800 and len({(r['language'],r['strategy'],r['nu'],r['a'],r['b']) for r in rows})==len(rows)
    summary = []
    for language in LANGUAGES:
        for nu in ORDERS:
            pair = [[r for r in rows if (r['language'],r['nu'],r['strategy'])==(language,nu,s)]
                    for s in STRATEGIES]
            assert len(pair[0])==len(pair[1])==1600
            assert [(r['a'],r['b']) for r in pair[0]]==[(r['a'],r['b']) for r in pair[1]]
            matched = [(r['microseconds'],s['microseconds']) for r,s in zip(*pair)
                       if r['setting'] is not None and s['setting'] is not None]
            summary.append(dict(language=language,nu=nu,paired_qualified=len(matched),
                                failed=[sum(r['setting'] is None for r in group) for group in pair],
                                ratio_of_means=sum(r for r,s in matched)/sum(s for r,s in matched),
                                methods=[dict(Counter(r['method'] for r in group)) for group in pair]))
    manifest = dict(summary=summary, power_before=power_before, power_after=power_after,
                    publication_timing_allowed=power_before['ac_connected'] and power_after['ac_connected'],
                    reference='Section 5.1 audited defining Bessel integral, target relative error 1e-10',
                    timing='Median of five batches; routing, bound, calibration, references and accuracy checks excluded. Identical workloads shared between strategies.',
                    diagonal='Both strategies switch points with a=b otherwise assigned to AF to SE.',
                    rr='Original 939 Fortran branch and its separate Python and C++ ports; native language seed evaluation timed.',
                    unique_workloads=len(jobs), failures=failures, cpp_compile_command=cpp_command,
                    inputs_sha256={name:sha(BENCHMARK/name) for name in ('calibration.jsonl','audit.json')})
    (OUT/'timings.json').write_text(json.dumps(rows,allow_nan=False),encoding='utf-8')
    with (OUT/'timings.csv').open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (OUT/'summary.json').write_text(json.dumps(manifest,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps(summary,indent=2), flush=True)


if __name__=='__main__':
    main()
