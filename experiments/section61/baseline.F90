! Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
! Paper: Efficient numerical computation of the generalized Marcum
! Q-function with applications (2026).
! Preprint: https://arxiv.org/abs/2611.29209
! Citation: CITATION.cff / CITATION.bib at the repository root.
! Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

module cdf_baseline
  use iso_c_binding, only: dp => c_double, c_int, c_int64_t
#ifdef GAMMACHI
  use GammaCHI, only: cdfgamNC
#endif
  implicit none
contains
  subroutine baseline_values(a, nu, count, b, values, flags) bind(C)
    real(dp), value :: a, nu
    integer(c_int64_t), value :: count
    real(dp), intent(in) :: b(count)
    real(dp), intent(out) :: values(count)
    integer(c_int), intent(out) :: flags(count)
    real(dp) :: q, df, lambda, threshold
#ifndef GAMMACHI
    real(dp) :: bound
    external :: cdfchn
#endif
    integer(c_int64_t) :: i
    df = 2.0_dp * nu
    lambda = a * a
    do i = 1, count
      if (b(i) == 0.0_dp) then
        values(i) = 0.0_dp
        flags(i) = 0
        cycle
      end if
      threshold = b(i) * b(i)
#ifdef GAMMACHI
      call cdfgamNC(2, df, lambda, threshold, values(i), q, flags(i))
#else
      call cdfchn(1, values(i), q, threshold, df, lambda, flags(i), bound)
#endif
    end do
  end subroutine

  function baseline_time(a, nu, count, b, values, flags, repeats) result(seconds) bind(C)
    real(dp), value :: a, nu
    integer(c_int64_t), value :: count, repeats
    real(dp), intent(in) :: b(count)
    real(dp), intent(out) :: values(count)
    integer(c_int), intent(out) :: flags(count)
    real(dp) :: seconds
    real(dp), volatile :: checksum
    integer(c_int64_t) :: i, start, finish, rate
    checksum = 0.0_dp
    call system_clock(start,rate)
    do i = 1, repeats
      call baseline_values(a, nu, count, b, values, flags)
      checksum = checksum + values(count)
    end do
    call system_clock(finish)
    seconds = real(finish-start,dp)/real(rate,dp)/repeats
  end function
end module
