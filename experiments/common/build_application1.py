# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Isolated Application 1 libraries. Never writes upstream or older libraries.

Original 939 qser/pser bodies retain directions and numerical handling.
Recorded generated edits expose fixed counts and memoize identical complete
log-Gamma / log(x) / exp(-x) values during one approved batch only.
"""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

from source939 import ROOT, UPSTREAM, verify_sources
from build_rebuilt import fortran_text, cpp_text


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generated_fortran(text):
    text, edits = fortran_text(text)
    for name in ('qser', 'pser'):
        first = text.index('    SUBROUTINE '+name+'(')
        last = text.index('    END SUBROUTINE '+name, first)+len('    END SUBROUTINE '+name)
        body = text[first:last]
        body = body.replace('    USE IncgamFI', '    USE application1_cache\n    USE IncgamFI', 1)
        for old, new in [('loggam(mu+1.0_r8)', 'cdf_lgmu(mu,0)'),
                         ('loggam(mu+n+1.0_r8)', 'cdf_lgmu(mu,n)'),
                         ('loggam(n+1.0_r8)', 'cdf_lgn(n)'),
                         ('log(x)', 'cdf_logx(x)'), ('exp(-x)', 'cdf_expx(x)')]:
            if old in body:
                edits.append(dict(scope=name, original=old, replacement=new, occurrences=body.count(old)))
                body = body.replace(old, new)
        text = text[:first]+body+text[last:]
    return text, edits


def generated_cpp(text):
    text, edits = cpp_text(text)
    # Only original helpers and SE are compiled: no legacy batch/SE-RR API.
    text = text[:text.index('double ignega(')]
    cache = '''struct App1Cache {
 double lx,ex;std::array<double,1003> gm{},gn{};std::array<bool,1003> mr{},nr{};
 explicit App1Cache(double x):lx(std::log(x)),ex(std::exp(-x)){}
};
'''
    text = text.replace('struct Context {', cache+'struct Context {', 1)
    text = text.replace('int64_t fixed=-1; };', 'int64_t fixed=-1; App1Cache*cache=nullptr; };', 1)
    helpers = '''double cdf_logx(const Context&c){return c.cache?c.cache->lx:std::log(c.x);}
double cdf_expx(const Context&c){return c.cache?c.cache->ex:ex(-c.x);}
double cdf_lgmu(const Context&c,int n){
 if(!c.cache)return std::lgamma(c.nu+n+1.);
 if(!c.cache->mr[n]){c.cache->gm[n]=std::lgamma(c.nu+n+1.);c.cache->mr[n]=true;}
 return c.cache->gm[n];
}
double cdf_lgn(const Context&c,int n){
 if(!c.cache)return std::lgamma(n+1.);
 if(!c.cache->nr[n]){c.cache->gn[n]=std::lgamma(n+1.);c.cache->nr[n]=true;}
 return c.cache->gn[n];
}
'''
    at = text.index('void se(')
    body = text[at:]
    for old, new in [('std::lgamma(c.nu+1)', 'cdf_lgmu(c,0)'),
                     ('std::lgamma(c.nu+n+1)', 'cdf_lgmu(c,n)'),
                     ('std::lgamma(n+1.)', 'cdf_lgn(c,n)'),
                     ('std::log(c.x)', 'cdf_logx(c)'), ('ex(-c.x)', 'cdf_expx(c)')]:
        if old in body:
            edits.append(dict(scope='se', original=old, replacement=new, occurrences=body.count(old)))
            body = body.replace(old, new)
    return text[:at]+helpers+body+'\n}\n#include "application1.inc"\n', edits


def build():
    verified = verify_sources()
    src = UPSTREAM/'939/Fortran90/Src'
    for language, compiler_name in [('fortran', 'gfortran'), ('cpp', 'g++')]:
        compiler = shutil.which(compiler_name)
        if not compiler:
            raise RuntimeError('Compiler missing: '+compiler_name)
        folder = ROOT/'build/application1'/language
        folder.mkdir(parents=True, exist_ok=True)
        flags = ['-O3', '-fno-fast-math', '-ffp-contract=off', '-fPIC']
        commands = []
        def run(args):
            command = [str(x) for x in args]
            commands.append(command)
            result = subprocess.run(command, cwd=str(folder), stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, universal_newlines=True)
            if result.returncode:
                raise RuntimeError(result.stdout)
        library = folder/('application1.dll' if os.name == 'nt' else 'application1.so')
        if language == 'fortran':
            original = src/'marcumq.f90'
            text, edits = generated_fortran(original.read_text(encoding='utf-8'))
            generated = folder/'marcumq_application1.f90'
            generated.write_text(text, encoding='utf-8')
            flags += ['-frecursive', '-ffree-line-length-none', '-I', src, '-J', folder]
            sources = [src/'cons.f90', src/'incgam.f90', ROOT/'src/fortran/application1_cache.f90',
                       generated, ROOT/'src/fortran/application1.f90', ROOT/'src/fortran/application1_timer.f90']
            links = ['-static-libgfortran', '-static-libgcc']
        else:
            original = ROOT/'src/cpp/algorithm939.cpp'
            text, edits = generated_cpp(original.read_text(encoding='utf-8'))
            generated = folder/'se_application1.cpp'
            generated.write_text(text, encoding='utf-8')
            flags += ['-std=c++17', '-I', ROOT/'src/cpp', '-I', ROOT/'src/cpp']
            sources = [generated, ROOT/'src/cpp/application1_timer.cpp']
            links = ['-static-libgcc', '-static-libstdc++']
        objects = []
        for source in sources:
            obj = folder/(source.stem+'.o')
            run([compiler]+flags+['-c', source, '-o', obj])
            objects.append(obj)
        run([compiler, '-shared']+objects+links+['-o', library])
        extra = [ROOT/'src/cpp/application1.inc'] if language == 'cpp' else []
        manifest = dict(language=language, commands=commands, changes=edits,
                        source_hashes={str(p.relative_to(ROOT)): sha(p) for p in sources+extra+[original]},
                        generated_sha256=sha(generated), library_sha256=sha(library),
                        upstream_verified=verified,
                        authorization='User approved Application 1 plan 2026-10-30: direct CDF, both batch optimized, shared initialization, external accuracy calibration.',
                        cache_scope='One distribution per batch; resets every call; scalar mode disables all shared caches',
                        numerical_changes='P output and b=0 boundary; no new recurrence, stopping, scaling, compensation, clipping, or fallback')
        (folder/'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
        print(language+': '+str(library), flush=True)
    verify_sources()


if __name__ == '__main__':
    build()
