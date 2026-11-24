# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Floating-point setup used by the published timing drivers."""
import ctypes
from marcum._paths import ROOT
def load():
    lib = ctypes.CDLL(str(ROOT/'build/fp_environment.dll'))
    lib.million_set_binary64 = lib.hybrid_set_binary64
    return lib, {}
