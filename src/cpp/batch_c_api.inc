// Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
// Paper: Efficient numerical computation of the generalized Marcum
// Q-function with applications (2026).
// Preprint: https://arxiv.org/abs/2611.29209
// Citation: CITATION.cff / CITATION.bib at the repository root.
// Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

extern "C" int components_batch_plan_size_c(){return int(sizeof(components_batch_plan));}
extern "C" void components_triples_c(int method,int64_t count,const double*a,const double*b,const double*nu,int target,double atol,double rtol,int64_t limit,components_result*out){
 for(int64_t i=0;i<count;i++)components_eval_c(method,a[i],b[i],nu[i],target,atol,rtol,limit,0,10,&out[i]);
}
extern "C" void components_batch_prepare_c(int method,double a,double nu,double atol,double rtol,int64_t limit,double bmin,double bmax,components_batch_plan*plan){
 try{prepare_batch(method,a,nu,atol,rtol,limit,bmin,bmax,*plan);}catch(const std::exception&){plan->status=5;}
}
extern "C" void components_batch_apply_c(const components_batch_plan*plan,int64_t count,const double*b,components_result*out){apply_batch(*plan,count,b,out);}
extern "C" void components_batch_time_c(const components_batch_plan*plan,int64_t count,const double*b,int64_t repeats,int mode,double*seconds){
 if(count<0||repeats<1){*seconds=0;return;}std::vector<components_result> out(static_cast<size_t>(count));components_batch_plan scratch{};volatile double checksum=0;
 auto start=std::chrono::steady_clock::now();
 for(int64_t j=0;j<repeats;j++){
  const components_batch_plan*current=plan;
  if(mode!=0){components_batch_prepare_c(int(plan->method),plan->a,plan->nu,plan->atol,plan->rtol,plan->max_terms,plan->bmin,plan->bmax,&scratch);current=&scratch;}
  if(mode!=2){apply_batch(*current,count,b,out.data());if(count&&std::isfinite(out.back().value))checksum+=out.back().value;}
  else checksum+=current->se_count+current->serr_count;
 }
 *seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
}
