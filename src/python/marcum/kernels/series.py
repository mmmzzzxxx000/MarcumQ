# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Source-939 helpers required by the fixed-setting ports."""
import math
def _starting_p(mu, x, y, epsilon, limit):
    mulnmu, lnx, lny, leps = mu*math.log(mu), math.log(x), math.log(y), math.log(epsilon)
    n = x+5 if x < 2 else 1.5*x
    def solve(n, mode):
        if n <= 0 or not math.isfinite(n):
            raise ArithmeticError("startingpser produced an unusable starting index")
        if mode == 0:
            return (n-leps)/(math.log(n)-lnx)
        if mode == 1:
            return (2*n-leps+mulnmu-mu*math.log(mu+n))/(math.log(n)-lnx-lny+math.log(mu+n))
        return (2*n-leps-y+mu*lny-mu*math.log(mu+n)+mu)/(math.log(n)-lnx-lny+math.log(mu+n))
    for mode in (0, 1 if y <= mu else 2):
        if mode and mu+n <= y:
            break
        previous = 0.
        for _ in range(200):
            if abs(n-previous) <= 1:
                break
            previous, n = n, solve(n, mode)
        else:
            raise ArithmeticError("startingpser iteration limit")
        if mode == 0:
            n += 1
    count = int(n)+1
    if count < 1:
        raise ArithmeticError("startingpser returned a nonpositive count")
    if count+1 > limit:
        raise OverflowError("startingpser needs more than max_terms")
    return count

