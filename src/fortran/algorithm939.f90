! Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
! Paper: Efficient numerical computation of the generalized Marcum
! Q-function with applications (2026).
! Preprint: https://arxiv.org/abs/2611.29209
! Citation: CITATION.cff / CITATION.bib at the repository root.
! Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

! Same fixed-count adapter with the ORIGINAL 939 RR component exposed.
! Recurrence steps/direction, continued fraction, fallback and seeds are intact.
module rebuild_adapter
 use, intrinsic :: iso_c_binding
 use MarcumQ, only: qser,pser,pqasyxy,MarcumPQtrap,prec,qrec,pqasymu,bench_setting
 use paper_algorithm1, only: paper_algorithm1_q
 implicit none
contains
 function rebuilt_q(method,a,b,nu,setting) result(q) bind(C,name='rebuilt_q')
  integer(c_int),value :: method
  real(c_double),value :: a,b,nu
  integer(c_int64_t),value :: setting
  real(c_double) :: p,q,x,y
  integer :: ierr
  if(method==3)then
   q=paper_algorithm1_q(a,b,nu,setting)
   return
  endif
  bench_setting=int(setting)
  x=a*a/2;y=b*b/2;ierr=0
  select case(method)
  case(0)
   if(y>x+nu)then
    call qser(nu,x,y,p,q,ierr)
   else
    call pser(nu,x,y,p,q,ierr)
   endif
  case(1)
   call pqasyxy(nu,x,y,p,q,ierr)
  case(2)
   call MarcumPQtrap(nu,x,y,p,q,ierr)
  case(5)
   call pqasymu(nu,x,y,p,q,ierr)
  case(4)
   if(y>x+nu)then
    call qrec(nu,x,y,p,q,ierr)
   else
    call prec(nu,x,y,p,q,ierr)
   endif
  end select
 end function
 function rebuilt_time(method,a,b,nu,setting,repeats) result(seconds) bind(C,name='rebuilt_time')
  integer(c_int),value :: method
  real(c_double),value :: a,b,nu
  integer(c_int64_t),value :: setting,repeats
  real(c_double) :: seconds
  real(c_double),volatile :: checksum
  integer(c_int64_t) :: k,t0,t1,rate
  checksum=0
  call system_clock(t0,rate)
  do k=1,repeats
   checksum=checksum+rebuilt_q(method,a,b,nu,setting)
  enddo
  call system_clock(t1)
  seconds=real(t1-t0,c_double)/real(rate,c_double)/repeats
 end function
end module
