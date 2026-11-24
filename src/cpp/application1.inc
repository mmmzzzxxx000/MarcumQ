// Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
// Paper: Efficient numerical computation of the generalized Marcum
// Q-function with applications (2026).
// Preprint: https://arxiv.org/abs/2611.29209
// Citation: CITATION.cff / CITATION.bib at the repository root.
// Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

// Included after the isolated source-flow-preserving 939 SE port.
namespace {
std::vector<double> app1_u(double a,double nu,int64_t N){
 std::vector<double> u(N+1);u[0]=1.;u[1]=(a*a-2.*nu)/(4.*nu);
 for(int64_t n=2;n<=N;n++){
  double Den=4.*n*(nu+n-1.);
  u[n]=((a*a-2.*(nu+2.*n-2.))/Den)*u[n-1]-(1./Den)*u[n-2];
 }
 return u;
}
double app1_sum(double b,double nu,int64_t N,const std::vector<double>&u){
 double Sum=0.,PowerTerm=std::pow(b,2.*nu);
 for(int64_t n=0;n<=N;n++){
  Sum=Sum+u[n]*PowerTerm/(2.*(nu+n));PowerTerm=PowerTerm*(b*b);
 }
 return Sum;
}
double app1_scalar(int method,double a,double b,double nu,int64_t N,App1Cache*cache){
 if(b==0.)return 0.;
 if(method==1){auto u=app1_u(a,nu,N);return std::exp(-a*a/2.)/(std::pow(2.,nu-1.)*std::tgamma(nu))*app1_sum(b,nu,N,u);}
 Context c{a,b,nu,a*a/2.,b*b/2.,0.,0.,1e-17,1e-15,10000,10,true};
 c.fixed=N;c.cache=cache;components_result r{};se(c,0,r);return r.p;
}
}
extern "C" void app1_values(int method,int mode,double a,double nu,int64_t m,const double*bs,const int64_t*settings,double*out){
 if(mode==0){for(int64_t i=0;i<m;i++)out[i]=app1_scalar(method,a,bs[i],nu,settings[i],nullptr);}
 else if(method==0){
  App1Cache cache(a*a/2.);
  for(int64_t i=0;i<m;i++)out[i]=app1_scalar(method,a,bs[i],nu,settings[i],&cache);
 }else{
  int64_t N=*std::max_element(settings,settings+m);auto u=app1_u(a,nu,N);
  double factor=std::exp(-a*a/2.)/(std::pow(2.,nu-1.)*std::tgamma(nu));
  for(int64_t i=0;i<m;i++)out[i]=factor*app1_sum(bs[i],nu,settings[i],u);
 }
}
extern "C" void app1_serr_scan(int mode,double a,double nu,int64_t m,const double*bs,int64_t last,double*out){
 auto u=app1_u(a,nu,last);double factor=std::exp(-a*a/2.)/(std::pow(2.,nu-1.)*std::tgamma(nu));
 for(int64_t i=0;i<m;i++){
  double Sum=0.,PowerTerm=std::pow(bs[i],2.*nu);
  for(int64_t n=0;n<=last;n++){
   Sum=Sum+u[n]*PowerTerm/(2.*(nu+n));PowerTerm=PowerTerm*(bs[i]*bs[i]);
   out[i*(last+1)+n]=bs[i]==0.?0.:factor*Sum;
  }
 }
}
