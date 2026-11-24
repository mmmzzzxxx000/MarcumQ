! Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
! Paper: Efficient numerical computation of the generalized Marcum
! Q-function with applications (2026).
! Preprint: https://arxiv.org/abs/2611.29209
! Citation: CITATION.cff / CITATION.bib at the repository root.
! Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

module cdf_hybrid_fixed
  use iso_c_binding, only: dp => c_double, c_int, c_int64_t
  use application1_kernels, only: scalar_cdf, make_u, sum_u, app1_cached, app1_cache_begin
  implicit none
contains
  subroutine fixed_values(mode, a, nu, count, b, methods, orders, values, flags) bind(C)
    integer(c_int), value :: mode
    real(dp), value :: a, nu
    integer(c_int64_t), value :: count
    real(dp), intent(in) :: b(count)
    integer(c_int), intent(in) :: methods(count)
    integer(c_int64_t), intent(in) :: orders(count)
    real(dp), intent(out) :: values(count)
    integer(c_int), intent(out) :: flags(count)
    real(dp) :: factor
    real(dp), allocatable :: u(:)
    integer(c_int64_t) :: i, nmax

    flags = 0
    app1_cached = .false.
    if (mode /= 0) then
      nmax = maxval(orders, mask=methods==1)
      if (nmax >= 1) then
        allocate(u(0:nmax))
        call make_u(a, nu, nmax, u)
        factor = exp(-a*a/2.0_dp)/(2.0_dp**(nu-1.0_dp)*gamma(nu))
      end if
      call app1_cache_begin(a*a/2.0_dp)
    end if
    do i = 1, count
      if (mode /= 0 .and. methods(i) == 1) then
        values(i) = factor*sum_u(b(i), nu, orders(i), u)
      else
        values(i) = scalar_cdf(methods(i), a, b(i), nu, orders(i))
      end if
    end do
    app1_cached = .false.
  end subroutine

  function fixed_time(mode, a, nu, count, b, methods, orders, values, flags, repeats) result(seconds) bind(C)
    integer(c_int), value :: mode
    real(dp), value :: a, nu
    integer(c_int64_t), value :: count, repeats
    real(dp), intent(in) :: b(count)
    integer(c_int), intent(in) :: methods(count)
    integer(c_int64_t), intent(in) :: orders(count)
    real(dp), intent(out) :: values(count)
    integer(c_int), intent(out) :: flags(count)
    real(dp) :: seconds
    real(dp), volatile :: checksum
    integer(c_int64_t) :: i, start, finish, rate

    checksum = 0.0_dp
    call system_clock(start, rate)
    do i = 1, repeats
      call fixed_values(mode, a, nu, count, b, methods, orders, values, flags)
      checksum = checksum + values(count)
    end do
    call system_clock(finish)
    seconds = real(finish-start, dp)/real(rate, dp)/repeats
  end function
end module
