# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

import ctypes as ct
import hashlib, json, math, random, shutil, statistics, subprocess, sys
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.special import gammaincc, gammaln
from marcum._paths import ROOT
from common import power_status, write_csv, write_json
OUT=Path(__file__).resolve().parent
APP=ROOT/'build/application1/fortran'
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def exact_error(value, reference):
    return abs(mp.mpf(float(value))-reference)/abs(reference) if reference else mp.mpf(abs(float(value)))
def hybrid_mask(a, nu, b):
    y = b*b/2
    expected = np.zeros(len(b), dtype=bool)
    test = (b > 0.0) & (b <= 7.0)
    z = y[test]
    log_m = math.log(2.0) - a*a/4.0 + z + nu*np.log(z) - gammaln(nu+1.0)
    use_rr = log_m < math.log(1.0-2e-3)
    other = ~use_rr
    use_rr[other] = 2*(1+np.exp(log_m[other])) <= 1000*gammaincc(nu,z[other])
    expected[test] = use_rr
    return expected

