! Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
! Paper: Efficient numerical computation of the generalized Marcum
! Q-function with applications (2026).
! Preprint: https://arxiv.org/abs/2611.29209
! Citation: CITATION.cff / CITATION.bib at the repository root.
! Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

! Thin C adapter only. The upstream Algorithm 939 sources remain byte-identical.
module algorithm939_c_interface
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use MarcumQ, only: marcum
  implicit none
contains
  function algorithm939_abi_version_c() result(version) bind(C, name="algorithm939_abi_version_c")
    integer(c_int) :: version
    version = 1_c_int
  end function algorithm939_abi_version_c

  subroutine algorithm939_eval_c(mu, x, y, p, q, ierr) bind(C, name="algorithm939_eval_c")
    real(c_double), value, intent(in) :: mu, x, y
    real(c_double), intent(out) :: p, q
    integer(c_int), intent(out) :: ierr
    integer :: native_error
    p = 0.0_c_double
    q = 0.0_c_double
    ierr = 2_c_int
    if (.not.ieee_is_finite(mu)) return
    if (.not.ieee_is_finite(x)) return
    if (.not.ieee_is_finite(y)) return
    ! Prevalidate: upstream marcum resets ierr after its domain checks.
    if (mu < 1.0_c_double .or. mu > 10000.0_c_double) return
    if (x < 0.0_c_double .or. x > 10000.0_c_double) return
    if (y < 0.0_c_double .or. y > 10000.0_c_double) return
    if (y == 0.0_c_double) then
      p = 0.0_c_double
      q = 1.0_c_double
      ierr = 0_c_int
      return
    end if
    call marcum(mu, x, y, p, q, native_error)
    ierr = int(native_error, c_int)
  end subroutine algorithm939_eval_c
end module algorithm939_c_interface
