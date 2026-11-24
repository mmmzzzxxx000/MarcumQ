// Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
// Paper: Efficient numerical computation of the generalized Marcum
// Q-function with applications (2026).
// Preprint: https://arxiv.org/abs/2611.29209
// Citation: CITATION.cff / CITATION.bib at the repository root.
// Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

// C++ binary64 component ports. SE/AF/NI formulas and flow follow ACM 939,
// Gil, Segura, Temme (2014), DOI 10.1145/2591004. Retain upstream attribution.
// Original Fortran archive is unchanged. Gamma and scaled-erfc primitives
// below are independent C++ implementations, not calls into Fortran/Python.
#include "components.h"
#include "zeta_coefficients.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <chrono>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>
namespace {
constexpr double PI=3.14159265358979323846,U=1.11022302462515654042e-16,LN2=.69314718055994530942;
double nanv(){return std::numeric_limits<double>::quiet_NaN();}
double ex(double t){return t>709.782712893384?INFINITY:std::exp(t);}
struct Context { double a,b,nu,x,y,atol,rtol,delta,eps;int64_t limit;int refinements;bool source; };
components_result empty(){components_result r{};r.value=r.p=r.q=r.log_value=r.truncation_bound=r.roundoff_proxy=r.condition=r.sum_abs=nanv();r.iterations=r.evaluations=-1;r.status=5;return r;}
void finish(components_result &r,double p,double q,int status,int target){
 r.p=p;r.q=q;r.value=target==0?p:q;r.status=status;
 if((status==0||status==6||status==7)&&(!std::isfinite(p)||!std::isfinite(q)))r.status=5;
 else if((status==0||status==6||status==7)&&(p<0||p>1||q<0||q>1))r.status=2;
 if(r.value>0)r.log_value=std::log(r.value);
}
std::pair<double,double> gamma_pair(double alpha,double z){
 if(z==0)return {0,1}; if(!(alpha>0&&z>=0))throw std::runtime_error("gamma domain");
 double lf=alpha*std::log(z)-z-std::lgamma(alpha);
 if(z<alpha+1){
  double term=1/alpha,sum=term;
  for(int n=1;n<100000;n++){term*=z/(alpha+n);sum+=term;if(std::abs(term)<=std::abs(sum)*2*U){double p=ex(lf)*sum;return {p,1-p};}}
 }else{
  const double tiny=1e-300;double b=z+1-alpha,c=1/tiny,d=1/b,h=d;
  for(int n=1;n<100000;n++){double an=-double(n)*(n-alpha);b+=2;d=an*d+b;if(std::abs(d)<tiny)d=tiny;c=b+an/c;if(std::abs(c)<tiny)c=tiny;d=1/d;double change=d*c;h*=change;if(std::abs(change-1)<=4*U){double q=ex(lf)*h;return {1-q,q};}}
 }
 throw std::runtime_error("gamma iteration limit");
}
double erfcx_positive(double z){
 if(z<26)return std::exp(z*z)*std::erfc(z);
 double term=1,sum=1;
 for(int k=1;k<100;k++){term*=-(2*k-1)/(2*z*z);sum+=term;if(std::abs(term)<2*U*std::abs(sum))break;}
 return sum/(std::sqrt(PI)*z);
}
int starting_p(const Context&c){
 double n=c.x<2?c.x+5:1.5*c.x,lnx=std::log(c.x),lny=std::log(c.y),le=std::log(c.eps),mul=c.nu*std::log(c.nu);
 for(int pass=0;pass<2;pass++){
  if(pass==1&&c.nu+n<=c.y)break;
  int mode=pass==0?0:(c.y>c.nu?2:1);double previous=0;int it=0;
  while(std::abs(n-previous)>1&&it++<200){previous=n;if(!(n>0&&std::isfinite(n)))throw std::runtime_error("series start");
   if(mode==0)n=(n-le)/(std::log(n)-lnx);
   else if(mode==1)n=(2*n-le+mul-c.nu*std::log(c.nu+n))/(std::log(n)-lnx-lny+std::log(c.nu+n));
   else n=(2*n-le-c.y+c.nu*lny-c.nu*std::log(c.nu+n)+c.nu)/(std::log(n)-lnx-lny+std::log(c.nu+n));
  }
  if(it>=200||!std::isfinite(n))throw std::runtime_error("series start limit");if(pass==0)n+=1;
 }
 if(n>double(c.limit)-2)throw std::overflow_error("series max terms");return int(n)+1;
}
void se(const Context&c,int target,components_result&r){
 double p=0,q=0;int status=0;int64_t count=0;
 if(c.y>c.x+c.nu){
  double q0=q=gamma_pair(c.nu,c.y).second;
  double lh=c.nu*std::log(c.y)-c.y-std::lgamma(c.nu+1);
  if(lh>std::log(std::numeric_limits<double>::min()*10)&&c.x<100){
   double h=ex(lh);int n=0;int64_t limit=std::min<int64_t>(1000,c.limit-1);
   while(q!=0&&q0/q>c.delta&&n<limit){q0=c.x*(q0+h)/(n+1.);h=c.x*c.y*h/((n+1.)*(c.nu+n+1));q+=q0;n++;}
   if(q!=0&&q0/q>c.delta)status=1;q*=ex(-c.x);count=n+1;
  }else{
   double sum=0;bool done=false;int64_t k=0;
   for(;k<std::min<int64_t>(10000,c.limit);){double t=ex(-c.x+k*std::log(c.x)-std::lgamma(k+1.))*gamma_pair(c.nu+k,c.y).second;sum+=t;k++;
    if((sum==0&&k>150)||(sum>0&&t/sum<c.delta&&k>10)){done=true;break;}}
   q=sum;count=k;if(!done)status=1;
  }p=1-q;
 }else{
  int n=starting_p(c)+1;count=n+1;
  double lh=-c.x-c.y+n*std::log(c.x)+(n+c.nu)*std::log(c.y)-std::lgamma(c.nu+n+1)-std::lgamma(n+1.);
  if(lh<std::log(std::numeric_limits<double>::min()*10)){
   double sum=0;for(int k=n;k>0;k--)sum+=ex(k*std::log(c.x)-std::lgamma(k+1.))*gamma_pair(c.nu+k,c.y).first;
   p=(sum+gamma_pair(c.nu,c.y).first)*ex(-c.x);
  }else{double h=ex(lh),p1=gamma_pair(c.nu+n,c.y).first*ex(-c.x+n*std::log(c.x)-std::lgamma(n+1.));
   while(n>0){h*=n*(c.nu+n)/(c.x*c.y);p1=n*p1/c.x+h;p+=p1;n--;}}
  q=1-p;
 }if(p==0||q==0)status=3;r.iterations=count;finish(r,p,q,status,target);
}
double ignega(int n,double x,double delta){
 double alpha=.5-n;
 if(x>1.5){double p=0,q=(x-1-alpha)*(x+1-alpha),rr=4*(x+1-alpha),s=1-alpha,ro=0,t=1,g=1;
  for(int k=0;k<10000;k++){if(t/g<=delta)return g*ex(alpha*std::log(x))/(x+1-alpha);p+=s;q+=rr;rr+=8;s+=2;double tau=p*(1+ro);ro=tau/(q-tau);t*=ro;g+=t;}
 }else{double t=1,s=1/alpha;int k=1;for(int it=0;it<10000;it++){if(std::abs(t/s)<=delta){double g=std::sqrt(PI);for(int j=1;j<=n;j++)g/=.5-j;return ex(x)*(g-ex(alpha*std::log(x))*s);}t*=-x/k;s+=t/(k+alpha);k++;}}
 throw std::runtime_error("ignega limit");
}
void af(const Context&c,int target,components_result&r){
 if(c.x==c.y){r.status=4;return;}
 double xi=c.a*c.b,rho=c.b/c.a,s=c.y>=c.x?1:-1,z=(c.y-c.x)*(c.y-c.x)/(c.x+c.y+xi),rm=ex(c.nu*std::log(rho));
 double pq=.5*rm*erfcx_positive(std::sqrt(z))/std::sqrt(rho),v=2*c.nu-1,rhom=v*(rho-1),rhop=2*(rho+1),cc=s*rm/std::sqrt(8*PI),an=std::sqrt(xi);
 std::vector<double> bn(1,1);int limit=int(std::min<int64_t>(100,c.limit-1));
 while(std::abs(bn.back())>c.delta&&int(bn.size())-1<limit){int n=int(bn.size());an*=(4*c.nu*c.nu-(2*n-1)*(2*n-1))/(8*n*xi);bn.push_back(an*(rhom-n*rhop)/(rho*(v+2*n)));}
 int n0=int(bn.size())-1;if(n0==0){r.status=1;return;}int nr=std::min(int(z)+1,n0);std::vector<double> phi(n0+1);
 phi[nr]=ex((nr-.5)*std::log(z))*ignega(nr,z,c.delta);
 for(int n=nr+1;n<=n0;n++)phi[n]=(-z*phi[n-1]+1)/(n-.5);
 for(int n=nr-1;n>=1;n--)phi[n]=(1-(n+.5)*phi[n+1])/z;
 for(int n=1;n<=n0;n++){cc=-cc;pq+=cc*bn[n]*phi[n];}pq*=ex(-z);
 int status=std::abs(bn.back())<=c.delta?0:1;
 double p=s==1?1-pq:pq,q=s==1?pq:1-pq;if(p==0||q==0)status=3;r.iterations=n0+1;finish(r,p,q,status,target);
}
double xminsin(double x){
 if(std::abs(x)>1)return 6*(x-std::sin(x))/(x*x*x);
 const double fk[]={1.95088260487819821294,-.0244124470324439564863,.000145741981563655,-.0000005073893903402518,1.1556455068443e-9,-1.85522118416e-12,2.2117315e-15,-2.035e-18,1.5e-21};
 double u0=0,u1=0,u2=0,tt=2*(2*x*x-1);for(int k=8;k>=0;k--){u2=u1;u1=u0;u0=tt*u1-u2+fk[k];}return (u0-u2)/2;
}
double zeta(double x,double y){
 double z=y-x-1;if(std::abs(z)<.05){auto ck=zeta_coefficients(x);double z2=z/((2*x+1)*(2*x+1)),sum=1,power=1;
  for(int k=1;k<=10;k++){power*=z2;double t=ck[k]*power;sum+=t;if(std::abs(t)<=1e-15)break;}return -z/std::sqrt(2*x+1)*sum;}
 double w=std::sqrt(1+4*x*y),half=x+y-w-std::log(2*y/(1+w));if(half<0)throw std::runtime_error("zeta sign");double zz=std::sqrt(2*half);return x+1<y?-zz:zz;
}
void ni(const Context&c,int target,components_result&r){
 double xs=c.x/c.nu,ys=c.y/c.nu,xis2=4*xs*ys,wxis=std::sqrt(1+xis2);if(ys==xs+1){r.status=4;return;}
 int64_t evaluations=0;
 auto trapsum=[&](double h,double d){double b0=3;
  auto f=[&](double theta){evaluations++;if(theta>b0)return 0.;if(std::abs(theta)<1e-10){double rt=(1+wxis)/(2*ys);return rt/(1-rt)*ex(-c.nu*wxis*theta*theta*.5);}
   double t2=theta*theta,st=std::sin(theta),ct=std::cos(theta),ts=theta/st,s2=st*st,wx=std::sqrt(ts*ts+xis2),xm=xminsin(theta),p=xm*t2*ts/6;
   double term1=(p*(ts+1)-t2-s2*xis2)/(ct*wx+wxis);p=p*(1+(ts+1)/(wx+wxis))/(1+wxis);double ep=c.nu*(term1-std::log1p(p));
   if(ep<=std::log(1e-16)){b0=std::min(theta,b0);return 0.;}double rt=(ts+wx)/(2*ys);p=(2*theta*std::pow(std::sin(theta/2),2)-xm*t2*theta/6)/(2*ys*s2);double dr=p*(1+ts/wx);
   return ex(ep)*(dr*st+(ct-rt)*rt)/(rt*(rt-2*ct)+1);};
  double sum=d==0?f(0)/2:0,theta=d==0?h:d,end=d==0?3-h/2:3;while(theta<end&&theta<b0){sum+=f(theta);theta+=h;}return sum*h;};
 double h=3./8,v=1,pq=trapsum(h,0),eps=c.source?1e-13:std::max(1e-15,c.delta);int nc=0;
 while((v>eps&&nc<c.refinements)||nc<=2){nc++;double q=trapsum(h,h/2);v=q==0?0:std::abs(pq/q-1);h/=2;pq=(pq+q)/2;}
 double zz=zeta(xs,ys),scaled=pq*ex(-c.nu*.5*zz*zz)/PI;double p=zz<0?1-scaled:-scaled,q=zz<0?scaled:1+scaled;
 int status=v<=eps?0:1;if(p==0||q==0)status=3;r.iterations=nc;r.evaluations=evaluations;finish(r,p,q,status,target);
}
double logadd(double a,double b){if(a==-INFINITY)return b;if(b==-INFINITY)return a;double h=std::max(a,b);return h+std::log1p(ex(std::min(a,b)-h));}
void serr(const Context&c,int target,components_result&r){
 double ly=std::log(c.y),lg=std::lgamma(c.nu+1),base=-c.x/2+c.nu*ly-lg;
 double lm=c.nu>=1?base+c.y+std::min(0.,std::log(c.nu)-ly+std::log(-std::expm1(-c.y))):base+c.y+std::log(2-ex(-c.y));
 auto ltol=[&](double lv){return logadd(c.atol?std::log(c.atol):-INFINITY,c.rtol?std::log(c.rtol)+lv:-INFINITY);};
 double sb=target?ltol(0)-LN2-std::log1p(c.rtol):(c.atol?std::log(c.atol):-INFINITY);
 if(lm<=sb){finish(r,0,1,7,target);r.iterations=0;r.truncation_bound=ex(lm);r.roundoff_proxy=0;return;}
 double l0=-c.x+c.nu*ly-lg,p2=0,p1=0,t=1,sum=0,corr=0,ab=0;int exponent=0;
 double factor=8*U*(1+std::abs(c.x)+std::abs(c.nu*ly)+std::abs(lg))+4*U;
 bool fast=l0>-600&&l0<600&&base>-600&&base<600&&c.y<=32;
 double pref=fast?ex(l0):0,enext=fast?ex(base)*c.y:0,budgetdiv=2*(1+c.rtol);
 auto materialize=[&](double value){if(value==0)return 0.;double ls=l0+exponent*LN2;if(ls>-700&&ls<700)return value*ex(ls);return std::copysign(ex(std::log(std::abs(value))+ls),value);};
 for(int64_t n=0;n<c.limit;n++){
  if(n==1)t=c.y*((c.x-c.nu)/(c.nu+1))*p1;
  else if(n>=2)t=(c.y/n)*((c.x-c.nu-2*n+2)/(c.nu+n))*p1-(c.y/n)*(c.y/(c.nu+n))*((c.nu+n-2)/(c.nu+n-1))*p2;
  if(!std::isfinite(t)){r.status=5;return;}
  if(std::max({std::abs(t),std::abs(p1),std::abs(p2),std::abs(sum),std::abs(corr),ab})>std::ldexp(1.,400)){
   t=std::ldexp(t,-512);p1=std::ldexp(p1,-512);p2=std::ldexp(p2,-512);sum=std::ldexp(sum,-512);corr=std::ldexp(corr,-512);ab=std::ldexp(ab,-512);exponent+=512;fast=false;}
  double ns=sum+t;corr+=std::abs(sum)>=std::abs(t)?(sum-ns)+t:(t-ns)+sum;sum=ns;ab+=std::abs(t);
  double s=sum+corr,p=fast?s*pref:materialize(s),q=1-p,value=target?q:p;
  if(!std::isfinite(p)){finish(r,p,q,2,target);r.iterations=n+1;r.sum_abs=materialize(ab);return;}
  if(fast&&enext>0){
   double tf=INFINITY;if(n+2>c.y)tf=enext*(c.nu>=1?c.nu/(c.nu+n+1):2)/(1-c.y/(n+2));
   enext*=c.y/(n+2);
   if(tf>(c.atol+c.rtol*std::abs(value))/budgetdiv&&n+1<c.limit){p2=p1;p1=t;continue;}
  }else if(fast)fast=false;
  if(fast&&value!=0){
   double tail=INFINITY;
   if(n+2>c.y)tail=base+(c.nu>=1?std::log(c.nu)-std::log(c.nu+n+1):LN2)+(n+1)*ly-std::lgamma(n+2.)-std::log1p(-c.y/(n+2));
   double tol=c.atol+c.rtol*std::abs(value);
   bool met=std::isfinite(tail)&&tol>0&&tail<=std::log(tol/budgetdiv);
   if(!met&&n+1<c.limit){p2=p1;p1=t;continue;}
   double sigma=ab*pref,proxy=factor*sigma+(target?U*std::abs(value):0);
   int status=met?0:1;
   if(met&&(value<0||value>1||s<0||proxy>tol/2))status=2;
   finish(r,p,q,status,target);r.iterations=n+1;r.truncation_bound=ex(tail);
   r.roundoff_proxy=proxy;r.sum_abs=sigma;
   r.condition=target?(q?(1+sigma)/std::abs(q):INFINITY):sigma/std::abs(p);
   return;
  }
  double lp=s?std::log(std::abs(s))+l0+exponent*LN2:-INFINITY;
  double lv=target?(q?std::log(std::abs(q)):-INFINITY):lp,lab=ab?std::log(ab)+l0+exponent*LN2:-INFINITY;
  double tail=INFINITY;if(n+2>c.y)tail=base+(c.nu>=1?std::log(c.nu)-std::log(c.nu+n+1):LN2)+(n+1)*ly-std::lgamma(n+2.)-std::log1p(-c.y/(n+2));
  bool met=std::isfinite(tail)&&std::isfinite(lv)&&tail<=ltol(lv)-LN2-std::log1p(c.rtol);int status=1;
  if(!met&&n+1<c.limit){p2=p1;p1=t;continue;}
  double proxylog=std::log(factor)+lab;if(target&&value)proxylog=logadd(proxylog,std::log(U)+lv);
  if(met){if(!std::isfinite(value)||value<0||value>1||s<0)status=2;else if(value==0)status=target?2:3;else status=proxylog<=ltol(lv)-LN2?0:2;}
  finish(r,p,q,status,target);r.log_value=lv;r.truncation_bound=ex(tail);r.roundoff_proxy=ex(proxylog);r.sum_abs=materialize(ab);r.condition=ex((target?logadd(0,lab):lab)-lv);r.iterations=n+1;
  if(met)return;p2=p1;p1=t;
 }
}
void baseline(const Context&c,int target,components_result&r){
 if(c.nu<1||c.nu>10000||c.x>10000||c.y>10000){r.status=4;return;}
 Context fixed=c;fixed.source=true;fixed.delta=1e-17;fixed.eps=1e-15;fixed.refinements=10;
 double xi=c.a*c.b,w=std::sqrt(4*c.x+2*c.nu),low=c.x+c.nu-w,high=c.x+c.nu+w;
 if(c.x<30)se(fixed,target,r);
 else if(c.nu*c.nu<2*xi&&xi>30)af(fixed,target,r);
 else if((c.nu>=135&&low<=c.y&&c.y<=high)||(c.nu<135&&low<=c.y&&c.y<=high)){
  // The RR / large-order branches are outside the configured low-order paper
  // suite. Never replace an unported branch by NI while calling it 939.
  r.status=4;
 }else ni(fixed,target,r);
}
#include "batch_components.inc"
}
extern "C" int components_abi_version_c(){return 1;}
extern "C" int components_result_size_c(){return int(sizeof(components_result));}
extern "C" void components_eval_c(int method,double a,double b,double nu,int target,double atol,double rtol,int64_t max_terms,int source_defaults,int max_refinements,components_result*out){
 *out=empty();try{
  if(!(std::isfinite(a)&&std::isfinite(b)&&std::isfinite(nu)&&a>=0&&b>=0&&nu>0&&atol>=0&&rtol>=0&&(atol+rtol)>0&&max_terms>0&&max_terms<=1000000))return;
  if(b==0){finish(*out,0,1,6,target);out->iterations=0;return;}
  Context c{a,b,nu,.5*a*a,.5*b*b,atol,rtol,source_defaults?1e-17:std::min(.01,std::max(1e-17,rtol>0?rtol*.1:1e-13)),source_defaults?1e-15:std::min(.01,std::max(1e-17,rtol>0?rtol*.1:1e-15)),max_terms,max_refinements,bool(source_defaults)};
  if(!std::isfinite(c.x)||!std::isfinite(c.y)||c.y==0)return;
  if((method<3||method==4)&&(nu>10000||c.x>10000||c.y>10000)){out->status=4;return;}
  if(a==0){auto pq=gamma_pair(nu,c.y);finish(*out,pq.first,pq.second,0,target);out->iterations=0;return;}
  if(method==0)se(c,target,*out);else if(method==1)af(c,target,*out);else if(method==2)ni(c,target,*out);else if(method==3)serr(c,target,*out);
  else if(method==4)baseline(c,target,*out);
  else if(method==5){
   if(b<=7){serr(c,target,*out);if(out->status==0||out->status==6||out->status==7)return;int first=out->status;
    components_result recovery=empty();se(c,target,recovery);recovery.reserved=256+(first<<16);
    if(recovery.status==0||recovery.status==6||recovery.status==7){recovery.reserved|=2048;*out=recovery;return;}
    *out=empty();baseline(c,target,*out);out->reserved=256+(first<<16);}
   else baseline(c,target,*out);
   if(out->status!=0&&out->status!=6&&out->status!=7){int flags=out->reserved;components_result recovery=empty();ni(c,target,recovery);recovery.reserved=flags|1024;*out=recovery;}
  }else out->status=4;
 }catch(const std::exception&){out->status=5;}
}
extern "C" void components_time_c(int method,double a,double b,double nu,int target,double atol,double rtol,int64_t max_terms,int source_defaults,int max_refinements,int64_t repeats,double *seconds){
 if(repeats<1){*seconds=0;return;}
 components_result result;volatile double checksum=0;
 decltype(&components_eval_c) volatile function=components_eval_c;
 auto begin=std::chrono::steady_clock::now();
 for(int64_t k=0;k<repeats;k++){function(method,a,b,nu,target,atol,rtol,max_terms,source_defaults,max_refinements,&result);if(std::isfinite(result.value))checksum=checksum+result.value;}
 auto end=std::chrono::steady_clock::now();*seconds=std::chrono::duration<double>(end-begin).count();
}
#include "batch_c_api.inc"
