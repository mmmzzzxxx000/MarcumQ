! Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
! Paper: Efficient numerical computation of the generalized Marcum
! Q-function with applications (2026).
! Preprint: https://arxiv.org/abs/2611.29209
! Citation: CITATION.cff / CITATION.bib at the repository root.
! Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

! Application 1: only approved direct CDF output and shared preparation.
module application1_kernels
 use,intrinsic::iso_c_binding
 use MarcumQ,only:qser,pser,bench_setting
 use application1_cache
 implicit none
contains
 subroutine make_u(a,nu,N,u)
  real(c_double),intent(in)::a,nu
  integer(c_int64_t),intent(in)::N
  real(c_double),intent(out)::u(0:N)
  real(c_double)::Den
  integer(c_int64_t)::k
  u(0)=1d0;u(1)=(a*a-2d0*nu)/(4d0*nu)
  do k=2,N
   Den=4d0*k*(nu+k-1d0)
   u(k)=((a*a-2d0*(nu+2d0*k-2d0))/Den)*u(k-1)-(1d0/Den)*u(k-2)
  enddo
 end subroutine
 function sum_u(b,nu,N,u) result(s)
  real(c_double),intent(in)::b,nu,u(0:*)
  integer(c_int64_t),intent(in)::N
  real(c_double)::s,power
  integer(c_int64_t)::k
  s=0d0;power=b**(2d0*nu)
  do k=0,N
   s=s+u(k)*power/(2d0*(nu+k))
   power=power*(b*b)
  enddo
 end function
 function scalar_cdf(method,a,b,nu,N) result(p)
  integer(c_int),intent(in)::method
  real(c_double),intent(in)::a,b,nu
  integer(c_int64_t),intent(in)::N
  real(c_double)::p,q,x,y
  real(c_double),allocatable::u(:)
  integer::ierr
  if(b==0d0)then
   p=0d0;return
  endif
  if(method==1)then
   allocate(u(0:N));call make_u(a,nu,N,u)
   p=exp(-a*a/2d0)/(2d0**(nu-1d0)*gamma(nu))*sum_u(b,nu,N,u)
   deallocate(u)
  else
   x=a*a/2d0;y=b*b/2d0;bench_setting=int(N)
   if(y>x+nu)then
    call qser(nu,x,y,p,q,ierr)
   else
    call pser(nu,x,y,p,q,ierr)
   endif
  endif
 end function
 subroutine app1_values(method,mode,a,nu,m,bs,settings,out) bind(C,name='app1_values')
  integer(c_int),value::method,mode
  real(c_double),value::a,nu
  integer(c_int64_t),value::m
  real(c_double),intent(in)::bs(m)
  integer(c_int64_t),intent(in)::settings(m)
  real(c_double),intent(out)::out(m)
  real(c_double),allocatable::u(:)
  real(c_double)::factor
  integer(c_int64_t)::i,N
  app1_cached=.false.
  if(mode==0)then
   do i=1,m
    out(i)=scalar_cdf(method,a,bs(i),nu,settings(i))
   enddo
  else if(method==0)then
   call app1_cache_begin(a*a/2d0)
   do i=1,m
    out(i)=scalar_cdf(method,a,bs(i),nu,settings(i))
   enddo
   app1_cached=.false.
  else
   N=maxval(settings);allocate(u(0:N));call make_u(a,nu,N,u)
   factor=exp(-a*a/2d0)/(2d0**(nu-1d0)*gamma(nu))
   do i=1,m
    out(i)=factor*sum_u(bs(i),nu,settings(i),u)
   enddo
   deallocate(u)
  endif
 end subroutine
 subroutine app1_serr_scan(mode,a,nu,m,bs,last,out) bind(C,name='app1_serr_scan')
  integer(c_int),value::mode
  real(c_double),value::a,nu
  integer(c_int64_t),value::m,last
  real(c_double),intent(in)::bs(m)
  real(c_double),intent(out)::out(0:last,m)
  real(c_double)::u(0:last),factor,s,power
  integer(c_int64_t)::i,n
  call make_u(a,nu,last,u)
  factor=exp(-a*a/2d0)/(2d0**(nu-1d0)*gamma(nu))
  do i=1,m
   s=0d0;power=bs(i)**(2d0*nu)
   do n=0,last
    s=s+u(n)*power/(2d0*(nu+n));power=power*(bs(i)*bs(i))
    out(n,i)=factor*s
   enddo
   if(bs(i)==0d0)out(:,i)=0d0
  enddo
 end subroutine
end module
