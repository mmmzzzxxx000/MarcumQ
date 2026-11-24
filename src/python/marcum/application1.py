# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Approved Application 1: literal CDF kernels, scalar and shared preparation.

SE follows algorithm939.se (939 qser/pser); only output P is exposed.
SE-RR follows paper_algorithm1; output is C*Sum, as approved 2026-10-30.
No stopping, clipping, compensation, scaling, error policy, or hybrid fallback.
"""
import ctypes
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time

import numpy as np
from scipy import special
from .serr import coefficients, sum_series

from marcum._paths import ROOT
METHODS = ('se', 'se_rr')
MODES = ('scalar', 'batch')
TINY_LOG = math.log(sys.float_info.min*10)


def se_cdf(a, b, nu, N):
    """Fixed-count port of 939, returning its existing P rather than 1-P."""
    if b == 0.0:
        return 0.0
    x, y = a*a/2, b*b/2
    if y > x+nu:
        q0 = q = float(special.gammaincc(nu, y))
        lh = nu*math.log(y)-y-float(special.gammaln(nu+1))
        if lh > TINY_LOG and x < 100:
            h = math.exp(lh)
            for n in range(N):
                q0 = x*(q0+h)/(n+1.)
                h = x*y*h/((n+1.)*(nu+n+1))
                q += q0
            return 1.0-math.exp(-x)*q
        total = 0.0
        for k in range(N+1):
            total += math.exp(-x+k*math.log(x)-float(special.gammaln(k+1)))*float(special.gammaincc(nu+k, y))
        return 1.0-total
    n = N+1
    lh = -x-y+n*math.log(x)+(n+nu)*math.log(y)-float(special.gammaln(nu+n+1))-float(special.gammaln(n+1))
    if lh < TINY_LOG:
        total = 0.0
        for k in range(n, 0, -1):
            total += math.exp(k*math.log(x)-float(special.gammaln(k+1)))*float(special.gammainc(nu+k, y))
        return (total+float(special.gammainc(nu, y)))*math.exp(-x)
    h = math.exp(lh)
    p1 = float(special.gammainc(nu+n, y))*math.exp(-x+n*math.log(x)-float(special.gammaln(n+1)))
    p = 0.0
    while n > 0:
        h = h*n*(nu+n)/(x*y)
        p1 = n*p1/x+h
        p += p1
        n -= 1
    return p


def serr_cdf(a, b, nu, N):
    if b == 0.0:
        return 0.0
    u = coefficients(a, nu, N)
    return math.exp(-a*a/2)/(2.0**(nu-1)*math.gamma(nu))*sum_series(b, nu, u)


def se_batch(a, bs, nu, settings):
    """Vectorized original qser/pser directions, with per-point fixed counts.

    Reuse x, log(x), exp(-x), and log-Gamma values at identical arguments.
    Incomplete Gamma remains b-dependent. Unusual underflow paths execute the
    corresponding original scalar path, never a different algorithm.
    """
    x = a*a/2
    lx, ex = math.log(x), math.exp(-x)
    y = bs*bs/2
    out = np.zeros(len(bs))
    qids = np.flatnonzero((bs > 0) & (y > x+nu))
    pids = np.flatnonzero((bs > 0) & (y <= x+nu))
    if len(qids):
        z = y[qids]
        lh = nu*np.log(z)-z-float(special.gammaln(nu+1))
        regular = (lh > TINY_LOG) & (x < 100)
        for i in qids[~regular]:
            out[i] = se_cdf(a, float(bs[i]), nu, int(settings[i]))
        ids = qids[regular]
        z = y[ids]
        q0 = special.gammaincc(nu, z)
        q = q0.copy()
        h = np.exp(lh[regular])
        ns = settings[ids]
        for n in range(int(ns.max()) if len(ns) else 0):
            active = np.flatnonzero(ns > n)
            q0[active] = x*(q0[active]+h[active])/(n+1.)
            h[active] = x*z[active]*h[active]/((n+1.)*(nu+n+1))
            q[active] = q[active]+q0[active]
        out[ids] = 1.0-ex*q
    if len(pids):
        z = y[pids]
        ns = settings[pids]+1
        # Compute each repeated complete log-Gamma once per distribution.
        unique, inverse = np.unique(ns, return_inverse=True)
        gm = special.gammaln(nu+unique+1.)[inverse]
        gn = special.gammaln(unique+1.)[inverse]
        lh = -x-z+ns*lx+(ns+nu)*np.log(z)-gm-gn
        regular = lh >= TINY_LOG
        for i in pids[~regular]:
            out[i] = se_cdf(a, float(bs[i]), nu, int(settings[i]))
        ids = pids[regular]
        z, ns = z[regular], ns[regular]
        h = np.exp(lh[regular])
        p1 = special.gammainc(nu+ns, z)*np.exp(-x+ns*lx-gn[regular])
        total = np.zeros(len(ids))
        for n in range(int(ns.max()) if len(ns) else 0, 0, -1):
            active = np.flatnonzero(ns >= n)
            h[active] = h[active]*n*(nu+n)/(x*z[active])
            p1[active] = n*p1[active]/x+h[active]
            total[active] = total[active]+p1[active]
        out[ids] = total
    return out


def serr_batch(a, bs, nu, settings):
    u = coefficients(a, nu, max(1, int(settings.max())) )
    factor = math.exp(-a*a/2)/(2.0**(nu-1)*math.gamma(nu))
    total = np.zeros(len(bs))
    power = np.power(bs, 2.0*nu)
    for n in range(len(u)):
        ids = np.flatnonzero(settings >= n)
        total[ids] = total[ids]+u[n]*power[ids]/(2.0*(nu+n))
        power[ids] = power[ids]*(bs[ids]*bs[ids])
    return factor*total


def python_values(method, mode, a, bs, nu, settings):
    if mode == 'batch':
        return (se_batch if method == 'se' else serr_batch)(a, bs, nu, settings)
    fn = se_cdf if method == 'se' else serr_cdf
    return np.asarray([fn(a, float(b), nu, int(n)) for b, n in zip(bs, settings)])


class Runtime:
    """Validated array inputs stay outside all timers. Native timers inside DLL."""
    def __init__(self, language):
        self.language = language
        if language == 'python':
            self.manifest = {'language': language}
            return
        folder = ROOT/'build/application1'/language
        self.manifest = json.loads((folder/'manifest.json').read_text(encoding='utf-8'))
        path = folder/('application1.dll' if os.name == 'nt' else 'application1.so')
        if hashlib.sha256(path.read_bytes()).hexdigest() != self.manifest['library_sha256']:
            raise RuntimeError('Application 1 library/manifest mismatch')
        self.dll = ctypes.CDLL(str(path))
        dp = ctypes.POINTER(ctypes.c_double)
        ip = ctypes.POINTER(ctypes.c_int64)
        args = [ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double,
                ctypes.c_int64, dp, ip, dp]
        self.dll.app1_values.argtypes = args
        self.dll.app1_values.restype = None
        self.dll.app1_time.argtypes = args+[ctypes.c_int64]
        self.dll.app1_time.restype = ctypes.c_double
        self.dll.app1_serr_scan.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double,
                                          ctypes.c_int64, dp, ctypes.c_int64, dp]
        self.dll.app1_serr_scan.restype = None

    def prepare(self, method, mode, a, bs, nu, settings):
        bs = np.ascontiguousarray(bs, dtype=np.float64)
        settings = np.ascontiguousarray(settings, dtype=np.int64)
        if bs.ndim != 1 or settings.shape != bs.shape or not len(bs):
            raise ValueError('Nonempty matching one-dimensional arrays required')
        if not (a > 0 and nu > 0 and np.all(np.isfinite(bs)) and np.all(bs >= 0)):
            raise ValueError('Application 1 requires a,nu>0 and finite b>=0')
        if np.any(settings < (1 if method == 'se_rr' else 0)) or np.any(settings > 1000):
            raise ValueError('Fixed counts outside declared Application 1 interval')
        output = np.empty(len(bs))
        if self.language == 'python':
            def values():
                return python_values(method, mode, a, bs, nu, settings)
            def timer(repeats):
                start = time.perf_counter()
                for _ in range(repeats):
                    result = python_values(method, mode, a, bs, nu, settings)
                seconds = (time.perf_counter()-start)/repeats
                return seconds
        else:
            args = (METHODS.index(method), MODES.index(mode), a, nu, len(bs),
                    bs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    settings.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)),
                    output.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
            def values():
                self.dll.app1_values(*args)
                return output.copy()
            def timer(repeats):
                return self.dll.app1_time(*(args+(repeats,)))
        return values, timer

    def serr_scan(self, mode, a, bs, nu, last):
        """External calibration only: exploit identical recurrence prefixes.

        Future coefficients never change a previous coefficient or partial sum.
        Returned columns are N=0..last (N=0 not a permitted timed kernel input).
        Every selected N is rechecked with the actual fixed-count kernel.
        """
        bs = np.ascontiguousarray(bs, dtype=np.float64)
        out = np.empty((len(bs), last+1))
        if self.language != 'python':
            self.dll.app1_serr_scan(MODES.index(mode), a, nu, len(bs),
                bs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), last,
                out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
            return out
        u = coefficients(a, nu, last)
        factor = math.exp(-a*a/2)/(2.0**(nu-1)*math.gamma(nu))
        with np.errstate(over='ignore', invalid='ignore', under='ignore'):
            if mode == 'batch':
                power = np.power(bs, 2.0*nu)
                total = np.zeros(len(bs))
                for n in range(last+1):
                    total = total+u[n]*power/(2.0*(nu+n))
                    power = power*(bs*bs)
                    out[:, n] = factor*total
            else:
                for i, b in enumerate(bs):
                    total = 0.0
                    power = float(b)**(2.0*nu)
                    for n in range(last+1):
                        total = total+u[n]*power/(2.0*(nu+n))
                        power = power*(float(b)*float(b))
                        out[i, n] = factor*total
        out[bs == 0, :] = 0.0
        return out
