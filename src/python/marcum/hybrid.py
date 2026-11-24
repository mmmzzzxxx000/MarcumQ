# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Fixed-work Python, C++ and Fortran hybrid entry points.

Method selection is outside numerical timing. The caller supplies the work
setting; the experimental drivers calibrate it against an independent reference.
"""
import ctypes
import importlib.util
from functools import lru_cache
import math
import sys
import numpy as np
from scipy.special import gammaincc, gammaln
from ._paths import ROOT

METHODS = ('SE', 'AF', 'NI', 'SE-RR', 'RR', 'LAE')


def select_method(a, b, nu, proposed=True):
    if not all(math.isfinite(v) and v > 0 for v in (a, b, nu)):
        raise ValueError('The manuscript hybrid requires finite a,b,nu > 0')
    x, y = a*a/2, b*b/2
    transition = x+nu-math.sqrt(4*x+2*nu) < y < x+nu+math.sqrt(4*x+2*nu)
    if proposed:
        lm = math.log(2)-x/2+y+nu*math.log(y)-gammaln(nu+1)
        if lm < math.log(.998):
            m = math.exp(lm)
            log_bound = math.log1p(m)-math.log1p(-m)
        else:
            tail = gammaincc(nu,y)
            log_bound = math.log(2)+float(np.logaddexp(0,lm))-math.log(tail) if tail > 0 else math.inf
        if b <= 7 and log_bound <= math.log(1000):
            return 'SE-RR'
        if (b > 7 and a <= 7) or b <= 7:
            return 'SE'
    elif x < 30:
        return 'SE'
    if a*b > 30 and nu*nu < 2*a*b:
        return 'SE' if a == b else 'AF'
    if transition:
        return 'RR' if nu < 135 else 'LAE'
    return 'NI'


@lru_cache(maxsize=1)
def _library():
    folder = ROOT/'experiments/common'
    sys.path.insert(0,str(ROOT/'experiments/common'))
    spec = importlib.util.spec_from_file_location('marcum_hybrid_build',folder/'build_hybrid.py')
    builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)
    builder.build()
    library = ctypes.CDLL(str(ROOT/'build/hybrid/measure.dll'))
    library.rebuilt_q.argtypes = [ctypes.c_int]+[ctypes.c_double]*3+[ctypes.c_int64]
    library.rebuilt_q.restype = ctypes.c_double
    library.marcum_fortran_hybrid_q.argtypes = [ctypes.c_double]*3+[ctypes.c_int64,ctypes.c_int]
    library.marcum_fortran_hybrid_q.restype = ctypes.c_double
    return library


@lru_cache(maxsize=1)
def _cpp_library():
    _library()
    library = ctypes.CDLL(str(ROOT/'build/hybrid/hybrid_cpp.dll'))
    library.marcum_cpp_hybrid_q.argtypes = [ctypes.c_double]*3+[ctypes.c_int64,ctypes.c_int]
    library.marcum_cpp_hybrid_q.restype = ctypes.c_double
    return library


def hybrid_q(a, b, nu, setting, proposed=True, language='fortran'):
    """Evaluate the selected method at the externally supplied fixed work count.

This function does not certify accuracy or introduce a new stopping rule.
"""
    method = select_method(a,b,nu,proposed)
    if isinstance(setting,bool) or int(setting) != setting or setting < 0:
        raise ValueError('setting must be a nonnegative integer')
    if method == 'SE-RR' and setting < 1:
        raise ValueError('SE-RR requires N >= 1')
    if language == 'fortran':
        return _library().marcum_fortran_hybrid_q(a,b,nu,int(setting),int(proposed))
    if language == 'cpp':
        return _cpp_library().marcum_cpp_hybrid_q(a,b,nu,int(setting),int(proposed))
    if language != 'python':
        raise ValueError('language must be python, cpp or fortran')
    if method == 'SE-RR':
        from .serr import evaluate
        return evaluate(a,b,nu,int(setting))
    if method == 'RR':
        from .rr import kernel
        return kernel(a,b,nu)[0](int(setting))
    if method == 'LAE':
        return _library().rebuilt_q(5,a,b,nu,int(setting))
    from . import algorithm939
    return getattr(algorithm939,method.lower())(a,b,nu,int(setting))
