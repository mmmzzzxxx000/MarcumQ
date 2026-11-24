# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Fixed-setting Marcum Q kernels used in the revised manuscript."""
from .serr import evaluate as se_rr
from .rebuilt import kernel
__version__ = "1.0.0rc1"
__all__ = ["se_rr", "kernel"]
