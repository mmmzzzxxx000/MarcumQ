# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

import ctypes, hashlib, json
from marcum._paths import ROOT
def load(variant='rr'):
    folder=ROOT/'build/hybrid_timing_fortran_rr'
    record=json.loads((folder/'manifest.json').read_text(encoding='utf-8'))
    path=folder/'hybrid_timing.dll'
    assert hashlib.sha256(path.read_bytes()).hexdigest()==record['library_sha256']
    lib=ctypes.CDLL(str(path));sig=[ctypes.c_int]+[ctypes.c_double]*3+[ctypes.c_int64]
    lib.rebuilt_q.argtypes=sig;lib.rebuilt_q.restype=ctypes.c_double
    lib.rebuilt_time.argtypes=sig+[ctypes.c_int64];lib.rebuilt_time.restype=ctypes.c_double
    return lib,record