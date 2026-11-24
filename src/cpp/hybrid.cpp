// Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
// Paper: Efficient numerical computation of the generalized Marcum
// Q-function with applications (2026).
// Preprint: https://arxiv.org/abs/2611.29209
// Citation: CITATION.cff / CITATION.bib at the repository root.
// Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

// Manuscript hybrid rules with the approved AF-diagonal treatment.
// SE/AF/NI/SE-RR/RR reuse the existing C++ kernels. LAE calls the official
// Fortran component through a C bridge; it is not a newly translated formula.
#include "hybrid.hpp"
#include <cmath>
#include <limits>

extern "C" double rebuilt_q(int,double,double,double,std::int64_t);
extern "C" double rr_port_q(double,double,double,std::int64_t);
extern "C" double hybrid_gamma_q(double,double);
extern "C" double hybrid_lae_q(double,double,double,std::int64_t);

namespace marcum {
int select_method(double a,double b,double nu,bool proposed) {
    if (!(std::isfinite(a)&&std::isfinite(b)&&std::isfinite(nu)&&a>0&&b>0&&nu>0)) return -1;
    const double x=a*a/2, y=b*b/2, w=std::sqrt(4*x+2*nu);
    if (proposed) {
        const double lm=std::log(2.)-x/2+y+nu*std::log(y)-std::lgamma(nu+1);
        double log_bound;
        if (lm<std::log(.998)) {
            const double m=std::exp(lm);
            log_bound=std::log1p(m)-std::log1p(-m);
        } else {
            const double q=hybrid_gamma_q(nu,y);
            const double log_sum=lm>0 ? lm+std::log1p(std::exp(-lm)) : std::log1p(std::exp(lm));
            log_bound=q>0 ? std::log(2.)+log_sum-std::log(q) : std::numeric_limits<double>::infinity();
        }
        if (b<=7 && log_bound<=std::log(1000.)) return SERR;
        if ((b>7 && a<=7) || b<=7) return SE;
    } else if (x<30) return SE;
    if (a*b>30 && nu*nu<2*a*b) return a==b ? SE : AF;
    if (x+nu-w<y && y<x+nu+w) return nu<135 ? RR : LAE;
    return NI;
}

double hybrid_q(double a,double b,double nu,std::int64_t setting,bool proposed) {
    const int method=select_method(a,b,nu,proposed);
    if (method<0 || setting<0 || (method==SERR && setting<1))
        return std::numeric_limits<double>::quiet_NaN();
    if (method==RR) return rr_port_q(a,b,nu,setting);
    if (method==LAE) return hybrid_lae_q(a,b,nu,setting);
    return rebuilt_q(method,a,b,nu,setting);
}
}

extern "C" int marcum_cpp_select_method(double a,double b,double nu,int proposed) {
    return marcum::select_method(a,b,nu,proposed!=0);
}
extern "C" double marcum_cpp_hybrid_q(double a,double b,double nu,std::int64_t setting,int proposed) {
    return marcum::hybrid_q(a,b,nu,setting,proposed!=0);
}
