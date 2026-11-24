! Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
! Paper: Efficient numerical computation of the generalized Marcum
! Q-function with applications (2026).
! Preprint: https://arxiv.org/abs/2611.29209
! Citation: CITATION.cff / CITATION.bib at the repository root.
! Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

module outage_939
  use iso_c_binding, only: dp => c_double, c_int, c_int64_t
  use MarcumQ, only: pser, MarcumPQtrap, bench_setting
  implicit none
contains
  subroutine outage_939_value(branch, x, y, setting, p, ierr) bind(C)
    integer(c_int), value :: branch, setting
    real(dp), value :: x, y
    real(dp), intent(out) :: p
    integer(c_int), intent(out) :: ierr
    real(dp) :: mu, q

    mu = 1.0_dp
    ierr = 0
    bench_setting = setting
    if (branch == 0) then
      call pser(mu, x, y, p, q, ierr)
    else
      call MarcumPQtrap(mu, x, y, p, q, ierr)
    end if
  end subroutine

  function outage_939_time(branch, x, y, setting, repeats) result(seconds) bind(C)
    integer(c_int), value :: branch, setting
    real(dp), value :: x, y
    integer(c_int64_t), value :: repeats
    real(dp) :: seconds, mu, p, q
    real(dp), volatile :: checksum
    integer(c_int) :: ierr
    integer(c_int64_t) :: i, start, finish, rate

    mu = 1.0_dp
    ierr = 0
    bench_setting = setting
    checksum = 0.0_dp
    if (branch == 0) then
      call system_clock(start, rate)
      do i = 1, repeats
        call pser(mu, x, y, p, q, ierr)
        checksum = checksum + p
      end do
      call system_clock(finish)
    else
      call system_clock(start, rate)
      do i = 1, repeats
        call MarcumPQtrap(mu, x, y, p, q, ierr)
        checksum = checksum + p
      end do
      call system_clock(finish)
    end if
    seconds = real(finish-start, dp)/real(rate, dp)/repeats
  end function
end module
