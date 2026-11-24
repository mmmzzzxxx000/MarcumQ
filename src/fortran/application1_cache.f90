! Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
! Paper: Efficient numerical computation of the generalized Marcum
! Q-function with applications (2026).
! Preprint: https://arxiv.org/abs/2611.29209
! Citation: CITATION.cff / CITATION.bib at the repository root.
! Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

! Approved batch reuse of b-independent initialization in original 939 SE.
module application1_cache
 use Someconstants,only:r8
 use IncgamFI,only:loggam
 implicit none
 logical :: app1_cached=.false.,mu_ready(0:1002),n_ready(0:1002)
 real(r8) :: saved_lx,saved_ex,mu_lg(0:1002),n_lg(0:1002)
contains
 subroutine app1_cache_begin(x)
  real(r8),intent(in)::x
  app1_cached=.true.;mu_ready=.false.;n_ready=.false.
  saved_lx=log(x);saved_ex=exp(-x)
 end subroutine
 function cdf_logx(x) result(v)
  real(r8),intent(in)::x
  real(r8)::v
  if(app1_cached)then
   v=saved_lx
  else
   v=log(x)
  endif
 end function
 function cdf_expx(x) result(v)
  real(r8),intent(in)::x
  real(r8)::v
  if(app1_cached)then
   v=saved_ex
  else
   v=exp(-x)
  endif
 end function
 function cdf_lgmu(mu,n) result(v)
  real(r8),intent(in)::mu
  integer,intent(in)::n
  real(r8)::v
  if(app1_cached)then
   if(.not.mu_ready(n))then
    mu_lg(n)=loggam(mu+n+1.0_r8);mu_ready(n)=.true.
   endif
   v=mu_lg(n)
  else
   v=loggam(mu+n+1.0_r8)
  endif
 end function
 function cdf_lgn(n) result(v)
  integer,intent(in)::n
  real(r8)::v
  if(app1_cached)then
   if(.not.n_ready(n))then
    n_lg(n)=loggam(n+1.0_r8);n_ready(n)=.true.
   endif
   v=n_lg(n)
  else
   v=loggam(n+1.0_r8)
  endif
 end function
end module
