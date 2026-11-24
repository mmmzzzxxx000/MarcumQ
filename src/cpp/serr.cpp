// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Zhong-Xuan Mao and Jing-Feng Tian
// Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
// Paper: Efficient numerical computation of the generalized Marcum
// Q-function with applications (2026).
// Preprint: https://arxiv.org/abs/2611.29209
// Citation: CITATION.cff / CITATION.bib at the repository root.
// Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

// Fixed-truncation Algorithm 1 with an equivalent formula for nu > 100.
#include <cmath>
#include <cstdint>
#include <vector>

namespace {
double evaluate_large_order(double a, double b, double nu, int64_t N) {
    const double a2 = a*a;
    const double b2 = b*b;

    // Generate b^(2n)*u_n directly, then sum in ascending order.
    std::vector<double> v(static_cast<std::size_t>(N+1));
    v[0] = 1.0;
    v[1] = b2*(a2 - 2.0*nu)/(4.0*nu);
    for (int64_t n = 2; n <= N; ++n) {
        const double denominator = 4.0*n*(nu + n - 1.0);
        v[n] = (b2*(a2 - 2.0*(nu + 2.0*n - 2.0))/denominator)*v[n-1]
               - (b2*b2/denominator)*v[n-2];
    }

    double total = 0.0;
    for (int64_t n = 0; n <= N; ++n) {
        total += (nu/(nu+n))*v[n];
    }

    const double y = b2/2.0;
    if (y == 0.0 || total == 0.0) return 1.0;
    const double log_factor = -a2/2.0 + nu*std::log(y) - std::lgamma(nu + 1.0);
    return 1.0 - std::copysign(std::exp(log_factor + std::log(std::abs(total))), total);
}
}

extern "C" double paper_algorithm1_q(double a, double b, double nu, int64_t N) {
    if (nu > 100.0) return evaluate_large_order(a, b, nu, N);
    std::vector<double> u(static_cast<std::size_t>(N+1));
    u[0] = 1.0;
    u[1] = (a*a-2.0*nu)/(4.0*nu);
    for (int64_t n=2; n<=N; ++n) {
        double Den = 4.0*n*(nu+n-1.0);
        u[n] = ((a*a-2.0*(nu+2.0*n-2.0))/Den)*u[n-1] - (1.0/Den)*u[n-2];
    }
    double Sum = 0.0;
    double PowerTerm = std::pow(b, 2.0*nu);
    for (int64_t n=0; n<=N; ++n) {
        Sum = Sum + u[n]*PowerTerm/(2.0*(nu+n));
        PowerTerm = PowerTerm*(b*b);
    }
    return 1.0-std::exp(-a*a/2.0)/(std::pow(2.0,nu-1.0)*std::tgamma(nu))*Sum;
}
