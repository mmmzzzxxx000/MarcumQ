# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src/python"))
from marcum import se_rr, kernel
a,b,nu,N=2.0,2.0,1.0,40
print('SE-RR Python:',se_rr(a,b,nu,N))
for language in ('python','cpp','fortran'):
    evaluate,_=kernel('se_rr',a,b,nu,language)
    print(language,evaluate(N))
