# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Build the archived Section 5.2 random-timing kernels from original sources."""
import json
import shutil
import subprocess
from pathlib import Path
from marcum._paths import ROOT
from build_rebuilt import fortran_text, digest, replace_once
from source939 import UPSTREAM, verify_sources


def build():
    sources = verify_sources()
    folder = ROOT / 'build/hybrid'
    folder.mkdir(parents=True, exist_ok=True)
    here = ROOT / 'experiments/section52/random'
    compilers = {name: shutil.which(name) for name in ('g++', 'gfortran')}
    if not all(compilers.values()):
        raise RuntimeError('g++ and gfortran must be on PATH')
    upstream = UPSTREAM / '939/Fortran90/Src'
    original = upstream / 'marcumq.f90'
    generated, changes = fortran_text(original.read_text(encoding='utf-8'))
    generated = replace_once(generated,
        'PUBLIC :: marcum,qser,pser,pqasyxy,MarcumPQtrap,bench_setting',
        'PUBLIC :: marcum,qser,pser,pqasyxy,MarcumPQtrap,bench_setting,prec,qrec,pqasymu', changes)
    generated = replace_once(generated,
        'DO WHILE ((abs(bk/s)> 1.e-30_r8).AND.(k <= 16))',
        'DO WHILE (k <= 16 .AND. k <= bench_setting)', changes)
    generated_path = folder / 'marcumq_generated.f90'
    generated_path.write_text(generated, encoding='utf-8')
    dependencies = [upstream/'cons.f90', upstream/'incgam.f90', generated_path,
                    ROOT/'src/fortran/serr.f90', ROOT/'src/fortran/algorithm939.f90',
                    ROOT/'src/fortran/hybrid.f90', here/'measure.f90']
    cpp_source, cpp_changes = __import__('build_rebuilt').cpp_text(
        (ROOT/'src/cpp/algorithm939.cpp').read_text(encoding='utf-8'))
    cpp_source += '\n#include "rr.inc"\nextern "C" double hybrid_gamma_q(double a,double x){return gamma_pair(a,x).second;}\n'
    cpp_generated = folder/'algorithm939_generated.cpp'
    cpp_generated.write_text(cpp_source,encoding='utf-8')
    all_sources = dependencies + [here/'reference.cpp', cpp_generated, ROOT/'src/cpp/algorithm939.cpp',
        ROOT/'src/cpp/rr.inc', ROOT/'src/cpp/rebuild_adapter.inc', ROOT/'src/cpp/zeta_coefficients.hpp',
        ROOT/'src/cpp/serr.cpp', ROOT/'src/cpp/hybrid.cpp', ROOT/'src/cpp/hybrid.hpp']
    fingerprints = {str(p.relative_to(ROOT)): digest(p) for p in all_sources}
    outputs = [folder/'reference.dll', folder/'measure.dll', folder/'hybrid_cpp.dll']
    manifest = folder/'manifest.json'
    if manifest.exists() and all(p.exists() for p in outputs):
        old = json.loads(manifest.read_text(encoding='utf-8'))
        if old['source_sha256'] == fingerprints and all(old['library_sha256'][p.name] == digest(p) for p in outputs):
            return old
    flags = ['-O3', '-fno-fast-math', '-ffp-contract=off']
    commands = [
        [compilers['g++']] + flags + ['-std=c++17', '-shared', '-static', '-pthread',
            str(here/'reference.cpp'), '-o', str(outputs[0])],
        [compilers['gfortran']] + flags + ['-frecursive', '-ffree-line-length-none',
            '-I', str(upstream), '-J', str(folder), '-shared', '-static-libgfortran',
            '-static-libgcc'] + [str(p) for p in dependencies] + ['-o', str(outputs[1])],
        [compilers['g++']] + flags + ['-std=c++17','-shared','-static-libgcc','-static-libstdc++',
            '-I',str(ROOT/'src/cpp'),'-I',str(ROOT/'src/cpp'),str(cpp_generated),
            str(ROOT/'src/cpp/serr.cpp'),str(ROOT/'src/cpp/hybrid.cpp'),str(outputs[1]),
            '-o',str(outputs[2])]
    ]
    for command in commands:
        subprocess.run(command, cwd=str(folder), check=True)
    data = dict(commands=commands, source_sha256=fingerprints,
                library_sha256={p.name: digest(p) for p in outputs},
                original_source_sha256=sources, substitutions=changes,cpp_substitutions=cpp_changes,
                reference='Independent x87 extended-precision Bessel quadrature with a near-one probability bound',
                fixed_setting_scope='Only work counts and component visibility are exposed; upstream files remain unchanged.')
    manifest.write_text(json.dumps(data, indent=2), encoding='utf-8')
    return data
