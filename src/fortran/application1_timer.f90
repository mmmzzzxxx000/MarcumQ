! Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
! Paper: Efficient numerical computation of the generalized Marcum
! Q-function with applications (2026).
! Preprint: https://arxiv.org/abs/2611.29209
! Citation: CITATION.cff / CITATION.bib at the repository root.
! Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

module application1_timer
 use,intrinsic::iso_c_binding
 use application1_kernels,only:app1_values
 implicit none
contains
 function app1_time(method,mode,a,nu,m,bs,settings,out,repeats) result(seconds) bind(C,name='app1_time')
  integer(c_int),value::method,mode
  real(c_double),value::a,nu
  integer(c_int64_t),value::m,repeats
  real(c_double),intent(in)::bs(m)
  integer(c_int64_t),intent(in)::settings(m)
  real(c_double),intent(out)::out(m)
  real(c_double)::seconds
  real(c_double),volatile::checksum
  integer(c_int64_t)::r,t0,t1,rate
  checksum=0d0
  call system_clock(t0,rate)
  do r=1,repeats
   call app1_values(method,mode,a,nu,m,bs,settings,out)
   checksum=checksum+out(m)
  enddo
  call system_clock(t1)
  seconds=real(t1-t0,c_double)/real(rate,c_double)/repeats
 end function
end module
