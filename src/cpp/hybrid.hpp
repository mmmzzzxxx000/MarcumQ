// Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
// Paper: Efficient numerical computation of the generalized Marcum
// Q-function with applications (2026).
// Preprint: https://arxiv.org/abs/2611.29209
// Citation: CITATION.cff / CITATION.bib at the repository root.
// Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

#ifndef MARCUM_HYBRID_HPP
#define MARCUM_HYBRID_HPP
#include <cstdint>

namespace marcum {
enum Method { SE=0, AF=1, NI=2, SERR=3, RR=4, LAE=5 };
// Positive finite a,b,nu. Work counts are supplied by the caller.
int select_method(double a, double b, double nu, bool proposed=true);
double hybrid_q(double a, double b, double nu, std::int64_t setting, bool proposed=true);
}

extern "C" int marcum_cpp_select_method(double,double,double,int);
extern "C" double marcum_cpp_hybrid_q(double,double,double,std::int64_t,int);
#endif
