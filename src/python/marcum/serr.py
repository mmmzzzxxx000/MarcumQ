# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Zhong-Xuan Mao and Jing-Feng Tian
# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Fixed-truncation Algorithm 1 with an equivalent formula for nu > 100.

N is the LAST INDEX (N+1 terms). No accuracy policy enters this kernel.
"""
import math


def coefficients(a, nu, N):
    """Algorithm 1, lines 1--5. Caller supplies the paper's domain and N>=1."""
    u = [0.0] * (N + 1)
    u[0] = 1.0
    u[1] = (a*a - 2.0*nu) / (4.0*nu)
    for n in range(2, N + 1):
        Den = 4.0*n*(nu + n - 1.0)
        u[n] = ((a*a - 2.0*(nu + 2.0*n - 2.0))/Den)*u[n-1] - (1.0/Den)*u[n-2]
    return u


def sum_series(b, nu, u):
    """Algorithm 1, lines 6--10, including the final power update."""
    Sum = 0.0
    PowerTerm = b**(2.0*nu)
    for n in range(len(u)):
        Sum = Sum + u[n]*PowerTerm/(2.0*(nu+n))
        PowerTerm = PowerTerm*(b*b)
    return Sum


def _evaluate_large_order(a, b, nu, N):
    a2 = a*a
    b2 = b*b

    # Generate b^(2n)*u_n directly, then sum in ascending order.
    v = [0.0] * (N + 1)
    v[0] = 1.0
    v[1] = b2*(a2 - 2.0*nu)/(4.0*nu)
    for n in range(2, N + 1):
        denominator = 4.0*n*(nu + n - 1.0)
        v[n] = (b2*(a2 - 2.0*(nu + 2.0*n - 2.0))/denominator)*v[n-1] - (b2*b2/denominator)*v[n-2]

    total = 0.0
    for n in range(N + 1):
        total += (nu/(nu+n))*v[n]

    y = b2/2.0
    if y == 0.0 or total == 0.0:
        return 1.0
    log_factor = -a2/2.0 + nu*math.log(y) - math.lgamma(nu + 1.0)
    return 1.0 - math.copysign(math.exp(log_factor + math.log(abs(total))), total)


def evaluate(a, b, nu, N):
    """Algorithm 1, Q output only; ordinary binary64 operations, no clipping.

    N>=1 is required by the explicit initialization u[1] in the manuscript.
    A nonfinite result or exception must be recorded by the external
    experiment. For nu > 100, the equivalent expression uses log Gamma.
    """
    if nu > 100.0:
        return _evaluate_large_order(a, b, nu, N)
    u = coefficients(a, nu, N)
    Sum = sum_series(b, nu, u)
    return 1.0 - math.exp(-a*a/2.0)/(2.0**(nu-1.0)*math.gamma(nu))*Sum
