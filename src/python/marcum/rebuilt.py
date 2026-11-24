# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Four rebuilt raw kernels, separate from all legacy diagnostic APIs."""
import ctypes
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import time

METHODS=('se','af','ni','se_rr')
from marcum._paths import ROOT


@lru_cache(maxsize=2)
def native(language):
    if language not in ('cpp','fortran'):raise ValueError('Unknown native language')
    folder=ROOT/'build/rebuilt'/language
    path=folder/('rebuilt.dll' if os.name=='nt' else 'rebuilt.so')
    record=json.loads((folder/'manifest.json').read_text(encoding='utf-8'))
    if hashlib.sha256(path.read_bytes()).hexdigest()!=record['library_sha256']:
        raise RuntimeError('Rebuilt library does not match manifest')
    dll=ctypes.CDLL(str(path))
    fn=dll.rebuilt_q
    fn.argtypes=[ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_int64]
    fn.restype=ctypes.c_double
    timer=dll.rebuilt_time
    timer.argtypes=fn.argtypes+[ctypes.c_int64];timer.restype=ctypes.c_double
    return dll,fn,timer


def kernel(method,a,b,nu,language):
    """Return fixed-setting evaluator and raw numerical timer for this point.

    Repetitions, accuracy search and acceptance occur in the experiment.
    Native timers exclude Python/ctypes overhead. Python timers call the raw
    numerical function without Evaluation objects or post-hoc diagnostics.
    """
    if method not in METHODS:raise ValueError('Unknown method')
    if language=='python':
        if method=='se_rr':
            from .serr import evaluate as function
        else:
            from . import algorithm939
            function=getattr(algorithm939,method)
        def evaluate(setting):return function(a,b,nu,setting)
        def timer(setting,repeats):
            begin=time.perf_counter()
            for _ in range(repeats):function(a,b,nu,setting)
            return (time.perf_counter()-begin)/repeats
        return evaluate,timer
    dll,fn,ntimer=native(language)
    method_id=METHODS.index(method)
    def evaluate(setting):return fn(method_id,a,b,nu,setting)
    def timer(setting,repeats):return ntimer(method_id,a,b,nu,setting,repeats)
    return evaluate,timer
