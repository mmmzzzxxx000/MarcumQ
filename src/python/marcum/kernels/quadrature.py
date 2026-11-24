# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Source-939 helpers required by the fixed-setting ports."""
import math
from .zeta_coefficients import coefficients
def zeta(x, y):
    z = y-x-1
    if abs(z) < .05:
        ck = coefficients(x)
        z2, total, power = z/(2*x+1)**2, 1., 1.
        for k in range(1,11):
            power *= z2
            term = ck[k]*power
            total += term
            if abs(term) <= 1e-15:
                break
        return -z/math.sqrt(2*x+1)*total
    w = math.sqrt(1+4*x*y)
    half = math.fsum((x, y, -w, -math.log(2*y/(1+w))))
    if half < 0:
        raise ArithmeticError("negative squared zeta outside the local expansion")
    value = math.sqrt(2*half)
    return -value if x+1 < y else value


def xminsin(x):
    if abs(x) > 1:
        return 6*(x-math.sin(x))/(x*x*x)
    fk = [1.95088260487819821294, -.0244124470324439564863,
          .00014574198156365500, -.0000005073893903402518,
          .0000000011556455068443, -.00000000000185522118416,
          .0000000000000022117315, -.000000000000000002035, 1.5e-21]
    u0 = u1 = u2 = 0.
    tt = 2*(2*x*x-1)
    for c in reversed(fk):
        u2, u1 = u1, u0
        u0 = tt*u1-u2+c
    return (u0-u2)/2

