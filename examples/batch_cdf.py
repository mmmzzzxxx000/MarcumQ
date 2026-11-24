# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src/python"))
import numpy as np
from marcum.application1 import Runtime
b=np.linspace(0,5,1000)
runtime=Runtime('fortran')
evaluate,_=runtime.prepare('se_rr','batch',2.0,b,1.0,np.full(len(b),40))
print(evaluate())
