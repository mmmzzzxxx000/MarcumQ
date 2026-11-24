# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

import shutil, subprocess
from marcum._paths import ROOT
from build_rebuilt import build as raw
from build_939 import build as official
from build_application1 import build as cdf
from build_hybrid_rr import build as rr
from build_hybrid import build as hybrid
if __name__=='__main__':
    for tool in ('g++','gfortran'):
        if not shutil.which(tool): raise RuntimeError('Install '+tool+' and add its directory to PATH.')
    raw();official();cdf();rr();hybrid()
    subprocess.run(['g++','-O3','-shared',str(ROOT/'src/cpp/hybrid_fpu.cpp'),'-o',str(ROOT/'build/fp_environment.dll')],check=True)
