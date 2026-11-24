// Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
// Paper: Efficient numerical computation of the generalized Marcum
// Q-function with applications (2026).
// Preprint: https://arxiv.org/abs/2611.29209
// Citation: CITATION.cff / CITATION.bib at the repository root.
// Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

// Adapter for generated, source-flow-preserving 939 component ports.
// Fixed settings are selected outside timing using the shared reference.
extern "C" double paper_algorithm1_q(double,double,double,int64_t);
extern "C" double rebuilt_q(int method,double a,double b,double nu,int64_t setting) {
 if(method==3)return paper_algorithm1_q(a,b,nu,setting);
 Context c{a,b,nu,a*a/2,b*b/2,0.,0.,1e-17,1e-15,10000,10,true};
 c.fixed=setting;
 components_result r{};
 if(method==0)se(c,1,r);
 else if(method==1)af(c,1,r);
 else if(method==2)ni(c,1,r);
 else return nanv();
 return r.q;
}
extern "C" double rebuilt_time(int method,double a,double b,double nu,int64_t setting,int64_t repeats) {
 double (*volatile fn)(int,double,double,double,int64_t)=rebuilt_q;
 volatile double checksum=0.;
 auto t0=std::chrono::steady_clock::now();
 for(int64_t k=0;k<repeats;k++)checksum=checksum+fn(method,a,b,nu,setting);
 auto t1=std::chrono::steady_clock::now();
 return std::chrono::duration<double>(t1-t0).count()/repeats;
}
