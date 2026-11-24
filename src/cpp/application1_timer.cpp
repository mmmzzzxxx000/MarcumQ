// Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
// Paper: Efficient numerical computation of the generalized Marcum
// Q-function with applications (2026).
// Preprint: https://arxiv.org/abs/2611.29209
// Citation: CITATION.cff / CITATION.bib at the repository root.
// Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

#include <chrono>
#include <cstdint>
extern "C" void app1_values(int,int,double,double,int64_t,const double*,const int64_t*,double*);
extern "C" double app1_time(int method,int mode,double a,double nu,int64_t m,const double*bs,const int64_t*settings,double*out,int64_t repeats){
 volatile double checksum=0.;
 auto start=std::chrono::steady_clock::now();
 for(int64_t r=0;r<repeats;r++){
  app1_values(method,mode,a,nu,m,bs,settings,out);checksum=checksum+out[m-1];
 }
 return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()/repeats;
}
