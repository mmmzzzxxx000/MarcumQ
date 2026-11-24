// Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
// Paper: Efficient numerical computation of the generalized Marcum
// Q-function with applications (2026).
// Preprint: https://arxiv.org/abs/2611.29209
// Citation: CITATION.cff / CITATION.bib at the repository root.
// Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

// Native fixed-(a,nu) CDF plans. Included after scalar helpers, before the
// anonymous namespace closes. No reference/oracle is called during timing.
Context batch_context(const components_batch_plan&p,double b){
 return Context{p.a,b,p.nu,.5*p.a*p.a,.5*b*b,p.atol,p.rtol,
  std::min(.01,std::max(1e-17,p.rtol>0?p.rtol*.1:1e-13)),
  std::min(.01,std::max(1e-17,p.rtol>0?p.rtol*.1:1e-15)),p.max_terms,10,false};
}
void prepare_batch(int method,double a,double nu,double atol,double rtol,int64_t limit,double bmin,double bmax,components_batch_plan&p){
 p={};p.a=a;p.nu=nu;p.atol=atol;p.rtol=rtol;p.max_terms=limit;p.bmin=bmin;p.bmax=bmax;p.method=method;p.status=4;
 if(!(std::isfinite(a)&&std::isfinite(nu)&&std::isfinite(bmin)&&std::isfinite(bmax)&&a>=0&&nu>0&&bmin>=0&&bmax>=bmin&&atol>=0&&rtol>=0&&atol+rtol>0&&limit>0&&limit<=1000000&&method>=0&&method<=8))return;
 p.status=0;if(method<6)return;
 if(method==7||method==8){
  p.serr_count=std::min<int64_t>(161,limit);p.coefficients[0]=1;
  double x=.5*a*a;
  if(p.serr_count>1)p.coefficients[1]=(x-nu)/(nu+1);
  for(int n=2;n<p.serr_count;n++)p.coefficients[n]=(x-nu-2*n+2)/(n*(nu+n))*p.coefficients[n-1]-(nu+n-2)/(n*(nu+n-1)*(nu+n))*p.coefficients[n-2];
 }
 if(method==6||method==8){
  Context c=batch_context(p,bmin>0?bmin:bmax);double x=c.x;
  if(x>200||nu>10000)return; // Same-method scalar recovery remains available.
  components_result r=empty();if(bmin>0)se(c,0,r);else r.value=1;
  double budget=(atol+rtol*std::abs(r.value))/(4*(1+rtol));
  if(!(budget>0&&std::isfinite(budget)))return;
  double w=ex(-x);p.weights[0]=w;
  int last=int(std::min<int64_t>(512,limit-1));
  for(int k=0;k<=last;k++){
   if(k+2>x){double next=w*x/(k+1),tail=next/(1-x/(k+2));if(tail<=budget){p.se_count=k+1;break;}}
   if(k<last){w*=x/(k+1);p.weights[k+1]=w;}
  }
 }
}
void cached_se(const components_batch_plan&p,double b,components_result&r){
 Context c=batch_context(p,b);if(b==0){finish(r,0,1,6,0);r.iterations=0;return;}
 auto scalar=[&](){se(c,0,r);r.reserved|=512;};
 if(p.se_count<=0){scalar();return;}
 int K=int(p.se_count)-1;double lh=(p.nu+K)*std::log(c.y)-c.y-std::lgamma(p.nu+K+1);
 if(lh<-600||lh>600){scalar();return;}
 double g=gamma_pair(p.nu+K+1,c.y).first,h=ex(lh),sum=0,corr=0;
 for(int k=K;k>=0;k--){g+=h;double term=p.weights[k]*g,ns=sum+term;corr+=std::abs(sum)>=std::abs(term)?(sum-ns)+term:(term-ns)+sum;sum=ns;if(k)h*=(p.nu+k)/c.y;}
 double value=sum+corr;if(!(std::isfinite(value)&&value>=0&&value<=1)){scalar();return;}
 finish(r,value,1-value,value==0?3:0,0);r.iterations=K+1;
}
void cached_serr(const components_batch_plan&p,double b,components_result&r){
 if(b==0){finish(r,0,1,6,0);r.iterations=0;return;}
 Context c=batch_context(p,b);
 if(b>8||p.serr_count<=0||p.nu>10000){serr(c,0,r);r.reserved|=512;return;}
 double ly=std::log(c.y),lg=std::lgamma(p.nu+1),base=-c.x/2+p.nu*ly-lg;
 double pref=ex(-c.x+p.nu*ly-lg),enext=ex(base)*c.y,power=1,sum=0,corr=0,ab=0,value=0,tail=INFINITY;
 if(pref==0||!std::isfinite(pref)||enext==0){serr(c,0,r);r.reserved|=512;return;}
 double factor=4*U+8*U*(1+c.x+std::abs(p.nu*ly)+std::abs(lg));int status=1,n=0;
 for(;n<p.serr_count;n++){
  double term=p.coefficients[n]*power,ns=sum+term;corr+=std::abs(sum)>=std::abs(term)?(sum-ns)+term:(term-ns)+sum;sum=ns;ab+=std::abs(term);value=pref*(sum+corr);
  tail=INFINITY;if(n+2>c.y)tail=enext*(p.nu>=1?p.nu/(p.nu+n+1):2)/(1-c.y/(n+2));
  double budget=(p.atol+p.rtol*std::abs(value))/(2*(1+p.rtol));
  if(std::isfinite(tail)&&tail>0&&tail<=budget){status=(value<0||value>1||factor*ab*pref>(p.atol+p.rtol*std::abs(value))/2)?2:0;break;}
  if(!std::isfinite(value)){status=5;break;}
  power*=c.y;enext*=c.y/(n+2);
 }
 finish(r,value,1-value,status,0);r.iterations=std::min<int64_t>(n+1,p.serr_count);r.sum_abs=ab*pref;r.roundoff_proxy=factor*r.sum_abs;r.truncation_bound=tail;r.condition=value?r.sum_abs/std::abs(value):INFINITY;
}
void apply_batch(const components_batch_plan&p,int64_t count,const double*b,components_result*out){
 for(int64_t i=0;i<count;i++){
  out[i]=empty();if(p.status!=0||!std::isfinite(b[i])||b[i]<0||(b[i]>0&&b[i]<p.bmin)||b[i]>p.bmax){out[i].status=4;continue;}
  if(p.method<6){components_eval_c(int(p.method),p.a,b[i],p.nu,0,p.atol,p.rtol,p.max_terms,0,10,&out[i]);continue;}
  try{
   if(p.method==6)cached_se(p,b[i],out[i]);
   else if(p.method==7)cached_serr(p,b[i],out[i]);
   else{
    if(b[i]<=7){cached_serr(p,b[i],out[i]);if(out[i].status==0||out[i].status==6||out[i].status==7)continue;int reason=out[i].status;out[i]=empty();cached_se(p,b[i],out[i]);out[i].reserved|=256+(reason<<16);}
    else cached_se(p,b[i],out[i]);
   }
  }catch(const std::exception&){out[i].status=5;}
 }
}
