# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""A fixed-work evaluation using the manuscript's revised hybrid wrapper."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src/python'))
from marcum.hybrid import hybrid_q, select_method

a,b,nu = 8.,2.,1.
print('Selected method:',select_method(a,b,nu))
for language in ('python','cpp','fortran'):
    print(language, 'Q at externally supplied N=80:',hybrid_q(a,b,nu,80,language=language))
