// Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
// Paper: Efficient numerical computation of the generalized Marcum
// Q-function with applications (2026).
// Preprint: https://arxiv.org/abs/2611.29209
// Citation: CITATION.cff / CITATION.bib at the repository root.
// Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <atomic>
#include <thread>

namespace {
using LD = long double;
static_assert(std::numeric_limits<LD>::digits == 64 &&
              std::numeric_limits<LD>::max_exponent == 16384,
              "This reference requires x87 80-bit floating point.");
constexpr LD PI = 3.141592653589793238462643383279502884L;
struct Rule{std::vector<LD>x,w;};
Rule legendre(int n){Rule r;r.x.resize(n);r.w.resize(n);
 for(int i=0;i<(n+1)/2;i++){
  LD z=std::cos(PI*(i+.75L)/(n+.5L)),pp=0;
  for(int it=0;it<100;it++){LD p=1,p0=0;
   for(int j=1;j<=n;j++){LD p1=((2*j-1)*z*p-(j-1)*p0)/j;p0=p;p=p1;}
   pp=n*(z*p-p0)/(z*z-1);LD delta=p/pp;z-=delta;
   if(std::abs(delta)<2*std::numeric_limits<LD>::epsilon())break;
  }
  LD w=2/((1-z*z)*pp*pp);r.x[i]=-z;r.x[n-i-1]=z;r.w[i]=r.w[n-i-1]=w;
 }return r;
}
const Rule& rule(int n){static Rule a=legendre(16),b=legendre(32),c=legendre(64);
 return n==16?a:n==32?b:c;}

struct Density {
    LD a, nu;
    // Expand I_(nu-1) as a positive series. Start at its largest term
    // and sum outward so that the intermediate terms are no larger than one.
    LD operator()(LD t) const {
        const LD z = a*t;
        const LD x = z*z/4;
        const LD order = nu-1;
        const int mode = static_cast<int>(2*x/(std::sqrt(order*order+4*x)+order));
        LD sum = 1;
        LD term = 1;
        for (int k=mode; k>0; --k) {
            term *= k*(nu+k-1)/x;
            sum += term;
            if (term <= 1e-24L*sum) break;
        }
        term = 1;
        for (int k=mode+1;; ++k) {
            term *= x/(k*(nu+k-1));
            sum += term;
            if (term <= 1e-24L*sum) break;
        }
        const LD log_peak = (mode ? mode*std::log(x) : 0)
                          -std::lgamma(nu+mode)-std::lgamma(1.L+mode);
        return (2*nu-1)*std::log(t)-(nu-1)*std::log(2.L)
              -(t-a)*(t-a)/2-z + log_peak + std::log(sum);
    }
};

LD integrate(const Density& density, const std::vector<LD>& points, LD scale, int nodes) {
    const auto& r = rule(nodes);
    LD total = 0;
    for (size_t i=1; i<points.size(); ++i) {
        const LD mid = (points[i]+points[i-1])/2;
        const LD half = (points[i]-points[i-1])/2;
        LD sum = 0;
        for (int k=0; k<nodes; ++k)
            sum += r.w[k]*std::exp(density(mid+half*r.x[k])-scale);
        total += half*sum;
    }
    return total;
}

LD log_probability_bound(LD a, LD t, LD nu) {
    const LD z = (nu+std::sqrt(nu*nu+a*a*t*t))/(t*t);
    const LD s = (1-z)/2;
    return -s*t*t-nu*std::log(z)+a*a*s/z;
}

std::vector<LD> partition(LD b, LD centre, LD end, bool finer) {
    std::vector<LD> points{b,end};
    auto add = [&](LD x) { if (b<x && x<end) points.push_back(x); };
    for (LD offset : {-16.L,-8.L,-4.L,-2.L,0.L,2.L,4.L,8.L,16.L})
        add(centre+offset);
    const LD scale = std::max(1.L,b-centre);
    for (LD offset : {0.5L,2.L,8.L,32.L}) add(b+offset/scale);
    if (finer) {
        for (LD offset : {-12.L,-6.L,-3.L,-1.L,1.L,3.L,6.L,12.L})
            add(centre+offset);
        for (LD offset : {0.125L,1.L,4.L,16.L}) add(b+offset/scale);
    }
    std::sort(points.begin(),points.end());
    points.erase(std::unique(points.begin(),points.end()),points.end());
    return points;
}

void reference(double av, double bv, double nv, double* out) {
    const unsigned short word=0x037f;
    asm volatile("fldcw %0"::"m"(word));
    const LD a=av,b=bv,nu=nv;
    const LD centre = std::sqrt(a*a+2*nu);
    if (b<centre) {
        const LD bound=log_probability_bound(a,b,nu);
        if (bound<-38) {
            out[0]=1; out[1]=0; out[2]=0;
            out[3]=static_cast<double>(std::exp(bound)); out[4]=2;
            return;
        }
    }
    const LD anchor = std::max(centre,b);
    const LD end = anchor+32;
    const Density density{a,nu};
    const LD scale = density(anchor);
    const auto cuts=partition(b,centre,end,false);
    LD low=integrate(density,cuts,scale,16);
    LD high=integrate(density,cuts,scale,32);
    if (std::abs(high-low)>2e-14L*high) {
        const auto finer=partition(b,centre,end,true);
        low=integrate(density,finer,scale,32);
        high=integrate(density,finer,scale,64);
    }
    const LD log_q = scale+std::log(high);
    const LD relative_error=2*std::abs(high-low)/high+5e-14L
                           +std::exp(log_probability_bound(a,end,nu)-log_q);
    const LD q=std::exp(log_q);
    out[0]=static_cast<double>(q);
    out[1]=static_cast<double>(q-LD(out[0]));
    out[2]=static_cast<double>(log_q);
    out[3]=static_cast<double>(relative_error);
    out[4]=std::isfinite(log_q) && relative_error<1e-12L ? 1 : 0;
}
}

extern "C" void reference_batch(int n, const double* xyz, int threads, double* out) {
    std::atomic<int> next{0};
    auto work = [&] {
        for (int start = next.fetch_add(32); start < n; start = next.fetch_add(32))
            for (int i=start; i<std::min(n,start+32); ++i)
                reference(xyz[3*i],xyz[3*i+1],xyz[3*i+2],out+5*i);
    };
    std::vector<std::thread> workers;
    for (int i=1; i<threads; ++i) workers.emplace_back(work);
    work();
    for (auto& thread : workers) thread.join();
}
