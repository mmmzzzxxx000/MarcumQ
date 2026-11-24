! Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
! Paper: Efficient numerical computation of the generalized Marcum
! Q-function with applications (2026).
! Preprint: https://arxiv.org/abs/2611.29209
! Citation: CITATION.cff / CITATION.bib at the repository root.
! Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

! Manuscript hybrid selection and fixed-work evaluation. No accuracy search.
module marcum_hybrid
 use, intrinsic :: iso_c_binding
 use, intrinsic :: ieee_arithmetic
 use IncgamFI, only: incgam
 use rebuild_adapter, only: rebuilt_q
 implicit none
 private
 public :: select_method, hybrid_q, hybrid_lae_q
contains
 function select_method(a,b,nu,proposed) result(method) bind(C,name='marcum_fortran_select_method')
  real(c_double),value :: a,b,nu
  integer(c_int),value :: proposed
  integer(c_int) :: method
  real(c_double) :: x,y,w,lm,m,p,q,log_bound,log_sum
  integer :: ierr
  method=-1
  if (.not.(ieee_is_finite(a).and.ieee_is_finite(b).and.ieee_is_finite(nu))) return
  if (a<=0.or.b<=0.or.nu<=0) return
  x=a*a/2; y=b*b/2; w=sqrt(4*x+2*nu)
  if (proposed/=0) then
   lm=log(2.0_c_double)-x/2+y+nu*log(y)-log_gamma(nu+1)
   if (lm<log(.998_c_double)) then
    m=exp(lm)
    log_bound=log(1+m)-log(1-m)
   else
    call incgam(nu,y,p,q,ierr)
    if (lm>0) then
     log_sum=lm+log(1+exp(-lm))
    else
     log_sum=log(1+exp(lm))
    endif
    if (q>0) then
     log_bound=log(2.0_c_double)+log_sum-log(q)
    else
     log_bound=ieee_value(0.0_c_double,ieee_positive_inf)
    endif
   endif
   if (b<=7.and.log_bound<=log(1000.0_c_double)) then
    method=3
    return
   endif
   if ((b>7.and.a<=7).or.b<=7) then
    method=0
    return
   endif
  else if (x<30) then
   method=0
   return
  endif
  if (a*b>30.and.nu*nu<2*a*b) then
   method=1
   if (a==b) method=0
  else if (x+nu-w<y.and.y<x+nu+w) then
   method=4
   if (nu>=135) method=5
  else
   method=2
  endif
 end function

 function hybrid_q(a,b,nu,setting,proposed) result(q) bind(C,name='marcum_fortran_hybrid_q')
  real(c_double),value :: a,b,nu
  integer(c_int64_t),value :: setting
  integer(c_int),value :: proposed
  real(c_double) :: q
  integer(c_int) :: method
  method=select_method(a,b,nu,proposed)
  if (method<0.or.setting<0.or.(method==3.and.setting<1)) then
   q=ieee_value(0.0_c_double,ieee_quiet_nan)
   return
  endif
  q=rebuilt_q(method,a,b,nu,setting)
 end function

 function hybrid_lae_q(a,b,nu,setting) result(q) bind(C,name='hybrid_lae_q')
  real(c_double),value :: a,b,nu
  integer(c_int64_t),value :: setting
  real(c_double) :: q
  q=rebuilt_q(5_c_int,a,b,nu,setting)
 end function
end module
