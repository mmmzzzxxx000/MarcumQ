// Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
// Paper: Efficient numerical computation of the generalized Marcum
// Q-function with applications (2026).
// Preprint: https://arxiv.org/abs/2611.29209
// Citation: CITATION.cff / CITATION.bib at the repository root.
// Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

// Port of Algorithm 939's fc, prec and qrec, using the existing SE and NI ports.
namespace {
double rr_fc(double order,double z) {
 const double tiny=10*std::numeric_limits<double>::min();
 double ratio=tiny,c=tiny,d=0,change=0,b=2*order/z;
 int m=0;
 while(std::abs(change-1)>1e-15) {
  d=b+d;
  if(std::abs(d)<tiny)d=tiny;
  c=b+1/c;
  if(std::abs(c)<tiny)c=tiny;
  d=1/d;
  change=c*d;
  ratio*=change;
  b=2*(order+ ++m)/z;
 }
 return ratio;
}

double rr_seed(Context c,double order,int target) {
 c.nu=order;
 components_result out=empty();
 ni(c,target,out);
 return out.value;
}

double rr_prec(const Context& c) {
 const double limit=c.y-c.x+1+std::sqrt(2*(c.x+c.y)+1);
 const int steps=int(limit)+2-int(c.nu);
 const double shifted=c.nu+steps;
 double coef=std::sqrt(c.y/c.x)*rr_fc(shifted,2*std::sqrt(c.x*c.y));
 double p1=rr_seed(c,shifted+1,0),p0=rr_seed(c,shifted,0);
 for(int n=0;n<steps;++n) {
  const double p=((1+coef)*p0-p1)/coef;
  p1=p0;
  p0=p;
  coef=c.y/(shifted-n-1+c.x*coef);
 }
 return 1-p0;
}

double rr_qrec(const Context& c) {
 const double limit=c.y-c.x+(1-std::sqrt(2*(c.x+c.y)+1));
 if(limit<5) {
  if(c.x>=200)return rr_prec(c);
  components_result out=empty();
  se(c,1,out);
  return out.q;
 }
 const int steps=int(c.nu)-(int(limit)-1);
 const double shifted=c.nu-steps;
 std::array<double,301> coef;
 coef[0]=std::sqrt(c.y/c.x)*rr_fc(c.nu,2*std::sqrt(c.x*c.y));
 for(int n=1;n<=steps;++n)coef[n]=c.y/(c.nu-n+c.x*coef[n-1]);
 double q0=rr_seed(c,shifted-1,1),q1=rr_seed(c,shifted,1);
 for(int n=steps;n>=1;--n) {
  const double q=(1+coef[n])*q1-coef[n]*q0;
  q0=q1;
  q1=q;
 }
 return q1;
}
}

extern "C" double rr_port_q(double a,double b,double nu,int64_t setting) {
 Context c{a,b,nu,a*a/2,b*b/2,0.,0.,1e-17,1e-15,10000,10,true};
 c.fixed=setting;
 return c.y>c.x+c.nu?rr_qrec(c):rr_prec(c);
}

extern "C" double rr_port_time(double a,double b,double nu,int64_t setting,int64_t repeats) {
 double (*volatile fn)(double,double,double,int64_t)=rr_port_q;
 volatile double checksum=0.;
 const auto begin=std::chrono::steady_clock::now();
 for(int64_t k=0;k<repeats;++k)checksum=checksum+fn(a,b,nu,setting);
 const auto end=std::chrono::steady_clock::now();
 return std::chrono::duration<double>(end-begin).count()/repeats;
}
