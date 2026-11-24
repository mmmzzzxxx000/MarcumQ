# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Retime Section 6.2 with the original 939 branches and a shared error target."""

import ctypes as ct
import csv
import math
import os
import random
import shutil
import statistics
import subprocess
import sys

import mpmath as mp
import numpy as np

from pathlib import Path
from marcum._paths import ROOT
from marcum.application1 import Runtime
from common import power_status, write_csv, write_json
OUT=Path(__file__).resolve().parent
def relerr(value, reference): return abs(mp.mpf(float(value))-reference)/abs(reference)
from marcum.backends.algorithm939 import native_pair

sys.path.insert(0, str(ROOT / 'experiments/common'))
from build_rebuilt import UPSTREAM, digest, fortran_text, verify_sources

OUTPUT = ROOT / 'results/section62'


def build():
    sources = verify_sources()
    compiler = shutil.which('gfortran')
    if compiler is None:
        raise RuntimeError('gfortran not found')
    upstream = UPSTREAM / '939/Fortran90/Src'
    generated = OUTPUT / 'marcumq_fixed_count.f90'
    contents, edits = fortran_text((upstream / 'marcumq.f90').read_text(encoding='utf-8'))
    generated.write_text(contents, encoding='utf-8')
    dll = OUTPUT / 'outage_939.dll'
    command = [compiler, '-O3', '-fno-fast-math', '-ffp-contract=off', '-frecursive',
               '-ffree-line-length-none', '-I', str(upstream), '-shared',
               '-static-libgfortran', '-static-libgcc', str(upstream / 'cons.f90'),
               str(upstream / 'incgam.f90'), str(generated), str(OUT / 'retime_939.f90'),
               '-o', str(dll)]
    subprocess.run(command, cwd=OUTPUT, check=True)
    directory = os.add_dll_directory(os.path.dirname(compiler)) if hasattr(os, 'add_dll_directory') else None
    lib = ct.CDLL(str(dll))
    lib.runtime_directory = directory
    args = [ct.c_int, ct.c_double, ct.c_double, ct.c_int]
    lib.outage_939_value.argtypes = args + [ct.POINTER(ct.c_double), ct.POINTER(ct.c_int)]
    lib.outage_939_value.restype = None
    lib.outage_939_time.argtypes = args + [ct.c_int64]
    lib.outage_939_time.restype = ct.c_double
    return lib, dict(command=command, upstream_files=sources, changes=edits,
                     library_sha256=digest(dll), generated_sha256=digest(generated))


def value(lib, branch, x, y, setting):
    p, flag = ct.c_double(), ct.c_int()
    lib.outage_939_value(branch, x, y, setting, ct.byref(p), ct.byref(flag))
    return p.value, flag.value


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    lib, build_info = build()
    runtime = Runtime('fortran')
    with (ROOT / 'data/section62/references.csv').open(newline='', encoding='utf-8-sig') as handle:
        cases = list(csv.DictReader(handle))
    assert len(cases) == 99
    expected = {(float(k), float(t)) for k in list(range(5,36))+[29.5,30.5] for t in (0.1,0.01,0.001)}
    assert {(float(c['K']),float(c['gamma_th'])) for c in cases} == expected
    rows, jobs = [], []
    with mp.workdps(60):
        target = mp.mpf('1e-10')
        for case in cases:
            a, b = float(case['a']), float(case['b'])
            x, y = a*a/2, b*b/2
            p_ref = mp.mpf(case['reference'])
            assert y < x+1-math.sqrt(4*x+2) and 2*math.sqrt(x*y) < 30
            branch = int(x >= 30)

            original = native_pair(a, b, 1.0)
            default_p, default_flag = value(lib, branch, x, y, -1)
            assert original.ierr == default_flag == 0 and original.p == default_p, case
            default_error = relerr(default_p, p_ref)
            assert default_error <= target, case
            for setting in range(11 if branch else 1001):
                p, flag = value(lib, branch, x, y, setting)
                error_939 = relerr(p, p_ref)
                if flag == 0 and error_939 <= target:
                    break
            else:
                raise RuntimeError(f'No Algorithm 939 setting attained the target: {case}')

            bs = np.array([b])
            for order in range(1, 41):
                rr_value, rr_timer = runtime.prepare('se_rr', 'scalar', a, bs, 1.0, [order])
                p_rr = rr_value()[0]
                error_rr = relerr(p_rr, p_ref)
                if error_rr <= target:
                    break
            else:
                raise RuntimeError(f'No SE-RR order attained the target: {case}')
            assert order == int(case['timed_N']), case

            row = dict(K=float(case['K']), gamma_th=float(case['gamma_th']), a=a, b=b,
                       reference=case['reference'], algorithm939_branch='NI' if branch else 'SE',
                       algorithm939_setting=setting, algorithm939_value=p,
                       algorithm939_rel_error=float(error_939),
                       original939_default_rel_error=float(default_error),
                       se_rr_order=order, se_rr_value=float(p_rr), se_rr_rel_error=float(error_rr))
            rows.append(row)
            for name, timer in (
                ('algorithm939', lambda n, args=(branch, x, y, setting): lib.outage_939_time(*args, n)),
                ('se_rr', rr_timer),
            ):
                warm = timer(256)
                count = max(256, min(250000, round(0.012/warm)))
                jobs.append((row, name, timer, count, []))

    before = power_status()
    rng = random.Random(20261111)
    for _ in range(9):
        rng.shuffle(jobs)
        for _, _, timer, count, samples in jobs:
            samples.append(1e6*timer(count))
    after = power_status()
    for row, method, _, count, samples in jobs:
        row[method+'_us'] = statistics.median(samples)
        row[method+'_samples'] = samples
        row[method+'_repeats'] = count
    for row in rows:
        row['speedup'] = row['algorithm939_us']/row['se_rr_us']
    write_csv(OUTPUT / 'outage.csv', [{k: v for k, v in row.items() if not k.endswith('_samples')}
                                     for row in rows])
    write_json(OUTPUT / 'summary.json', dict(rows=rows, algorithm939_build=build_info,
               se_rr_build=runtime.manifest, power_before=before, power_after=after,
               target='1e-10', reference_source=str(ROOT / 'data/section62/references.csv'),
               timing='Nine shuffled native Fortran batches after external branch and setting selection.'))


if __name__ == '__main__':
    main()
