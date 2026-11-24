! Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
! Paper: Efficient numerical computation of the generalized Marcum
! Q-function with applications (2026).
! Preprint: https://arxiv.org/abs/2611.29209
! Citation: CITATION.cff / CITATION.bib at the repository root.
! Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

module cdf_baseline_fixed
  use iso_c_binding, only: dp => c_double, c_int, c_int64_t
  use GammaCHI, only: pser, qser, bench_setting, app1_cached, app1_cache_begin
  implicit none
contains
  subroutine fixed_values(mode, a, nu, count, b, tails, orders, values, flags) bind(C)
    integer(c_int), value :: mode
    real(dp), value :: a, nu
    integer(c_int64_t), value :: count
    real(dp), intent(in) :: b(count)
    integer(c_int), intent(in) :: tails(count)
    integer(c_int64_t), intent(in) :: orders(count)
    real(dp), intent(out) :: values(count)
    integer(c_int), intent(out) :: flags(count)
    real(dp) :: x, y, q
    integer(c_int64_t) :: i

    x = a*a/2.0_dp
    app1_cached = .false.
    if (mode /= 0) call app1_cache_begin(x)
    do i = 1, count
      flags(i) = 0
      if (b(i) == 0.0_dp) then
        values(i) = 0.0_dp
        cycle
      end if
      y = b(i)*b(i)/2.0_dp
      bench_setting = int(orders(i))
      if (tails(i) == 0) then
        call pser(nu, x, y, values(i), q, flags(i))
      else
        call qser(nu, x, y, values(i), q, flags(i))
      end if
    end do
    app1_cached = .false.
  end subroutine

  function fixed_time(mode, a, nu, count, b, tails, orders, values, flags, repeats) result(seconds) bind(C)
    integer(c_int), value :: mode
    real(dp), value :: a, nu
    integer(c_int64_t), value :: count, repeats
    real(dp), intent(in) :: b(count)
    integer(c_int), intent(in) :: tails(count)
    integer(c_int64_t), intent(in) :: orders(count)
    real(dp), intent(out) :: values(count)
    integer(c_int), intent(out) :: flags(count)
    real(dp) :: seconds
    real(dp), volatile :: checksum
    integer(c_int64_t) :: i, start, finish, rate

    checksum = 0.0_dp
    call system_clock(start, rate)
    do i = 1, repeats
      call fixed_values(mode, a, nu, count, b, tails, orders, values, flags)
      checksum = checksum + values(count)
    end do
    call system_clock(finish)
    seconds = real(finish-start, dp)/real(rate, dp)/repeats
  end function
end module
