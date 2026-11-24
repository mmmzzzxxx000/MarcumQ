! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 Zhong-Xuan Mao and Jing-Feng Tian
! Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
! Paper: Efficient numerical computation of the generalized Marcum
! Q-function with applications (2026).
! Preprint: https://arxiv.org/abs/2611.29209
! Citation: CITATION.cff / CITATION.bib at the repository root.
! Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

! Fixed-truncation Algorithm 1 with an equivalent formula for nu > 100.
module paper_algorithm1
 use, intrinsic :: iso_c_binding
 implicit none
 private :: evaluate_large_order
contains
 function paper_algorithm1_q(a,b,nu,N) result(Q) bind(C,name='paper_algorithm1_q')
  real(c_double),value :: a,b,nu
  integer(c_int64_t),value :: N
  real(c_double) :: Q,Den,Sum,PowerTerm
  real(c_double),allocatable :: u(:)
  integer(c_int64_t) :: k
  if (nu > 100.0_c_double) then
   Q=evaluate_large_order(a,b,nu,N)
   return
  end if
  allocate(u(0:N))
  u(0)=1.0_c_double
  u(1)=(a*a-2.0_c_double*nu)/(4.0_c_double*nu)
  do k=2,N
   Den=4.0_c_double*k*(nu+k-1.0_c_double)
   u(k)=((a*a-2.0_c_double*(nu+2.0_c_double*k-2.0_c_double))/Den)*u(k-1) &
         -(1.0_c_double/Den)*u(k-2)
  end do
  Sum=0.0_c_double
  PowerTerm=b**(2.0_c_double*nu)
  do k=0,N
   Sum=Sum+u(k)*PowerTerm/(2.0_c_double*(nu+k))
   PowerTerm=PowerTerm*(b*b)
  end do
  Q=1.0_c_double-exp(-a*a/2.0_c_double)/(2.0_c_double**(nu-1.0_c_double)*gamma(nu))*Sum
  deallocate(u)
 end function

 function evaluate_large_order(a,b,nu,N) result(Q)
  real(c_double),intent(in) :: a,b,nu
  integer(c_int64_t),intent(in) :: N
  real(c_double) :: Q,a2,b2,y,denominator,total,log_factor
  real(c_double),allocatable :: v(:)
  integer(c_int64_t) :: k

  a2=a*a
  b2=b*b

  ! Generate b^(2n)*u_n directly, then sum in ascending order.
  allocate(v(0:N))
  v(0)=1.0_c_double
  v(1)=b2*(a2 - 2.0_c_double*nu)/(4.0_c_double*nu)
  do k=2,N
   denominator=4.0_c_double*k*(nu + k - 1.0_c_double)
   v(k)=(b2*(a2 - 2.0_c_double*(nu + 2.0_c_double*k - 2.0_c_double))/denominator)*v(k-1) &
        -(b2*b2/denominator)*v(k-2)
  end do

  total=0.0_c_double
  do k=0,N
   total=total+(nu/(nu+k))*v(k)
  end do

  y=b2/2.0_c_double
  if (y == 0.0_c_double .or. total == 0.0_c_double) then
   Q=1.0_c_double
  else
   log_factor=-a2/2.0_c_double + nu*log(y) - log_gamma(nu + 1.0_c_double)
   Q=1.0_c_double-sign(exp(log_factor + log(abs(total))),total)
  end if
 end function
end module
