# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Source-939 helpers required by the fixed-setting ports."""
import math
def _ignega(n, x, delta):
    alpha = .5-n
    if x > 1.5:
        p, q = 0., (x-1-alpha)*(x+1-alpha)
        r, s, ro, t, g = 4*(x+1-alpha), 1-alpha, 0., 1., 1.
        for _ in range(10000):
            if t/g <= delta:
                return g*math.exp(alpha*math.log(x))/(x+1-alpha)
            p += s
            q += r
            r += 8
            s += 2
            tau = p*(1+ro)
            ro = tau/(q-tau)
            t *= ro
            g += t
    else:
        t, s, k = 1., 1/alpha, 1
        for _ in range(10000):
            if abs(t/s) <= delta:
                g = math.sqrt(math.pi)
                for j in range(1, n+1):
                    g /= .5-j
                return math.exp(x)*(g-math.exp(alpha*math.log(x))*s)
            t *= -x/k
            s += t/(k+alpha)
            k += 1
    raise ArithmeticError("ignega iteration limit")

