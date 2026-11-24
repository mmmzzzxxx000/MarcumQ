! Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
! Paper: Efficient numerical computation of the generalized Marcum
! Q-function with applications (2026).
! Preprint: https://arxiv.org/abs/2611.29209
! Citation: CITATION.cff / CITATION.bib at the repository root.
! Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

module random_hybrid_timing
    use iso_c_binding
    use ieee_arithmetic, only: ieee_is_finite
    use rebuild_adapter, only: rebuilt_q
    implicit none
    integer, parameter :: qp = selected_real_kind(30)
    real(c_double), volatile :: sink
contains
    function choose_setting(method,a,b,nu,hi,lo,bound) result(setting)
        integer(c_int), intent(in) :: method
        real(c_double), intent(in) :: a,b,nu,hi,lo,bound
        integer(c_int) :: setting,first,last,k
        real(c_double) :: value,x,y
        real(qp) :: ref,allowed
        ref = real(hi,qp)+real(lo,qp)
        allowed = 1.0e-10_qp*(ref-real(bound,qp))-real(bound,qp)
        setting=-1
        if (allowed <= 0) return
        first=0
        last=1000
        select case(method)
        case(1)
            first=1
            last=100
        case(2)
            last=10
        case(3)
            first=1
        case(5)
            last=16
        case(4)
            x=a*a/2
            y=b*b/2
            last=10
            if (y>x+nu .and. y-x+1-sqrt(2*(x+y)+1)<5 .and. x<200) last=1000
        end select
        do k=first,last
            value=rebuilt_q(method,a,b,nu,int(k,c_int64_t))
            if (.not.ieee_is_finite(value)) cycle
            if (abs(real(value,qp)-ref) <= allowed) then
                setting=k
                return
            end if
        end do
    end function

    subroutine calibrate(n,a,b,nu,hi,lo,err,old,new,old_n,new_n) bind(C)
        integer(c_int), value :: n
        real(c_double), intent(in) :: a(n),b(n),nu(n),hi(n),lo(n),err(n)
        integer(c_int), intent(in) :: old(n),new(n)
        integer(c_int), intent(out) :: old_n(n),new_n(n)
        integer :: i
        do i=1,n
            old_n(i)=choose_setting(old(i),a(i),b(i),nu(i),hi(i),lo(i),err(i))
            if (old(i)==new(i)) then
                new_n(i)=old_n(i)
            else
                new_n(i)=choose_setting(new(i),a(i),b(i),nu(i),hi(i),lo(i),err(i))
            end if
        end do
    end subroutine

    function batch_seconds(n,repeats,a,b,nu,method,setting) result(seconds) bind(C)
        integer(c_int), value :: n,repeats
        real(c_double), intent(in) :: a(n),b(n),nu(n)
        integer(c_int), intent(in) :: method(n),setting(n)
        real(c_double) :: seconds,checksum
        integer(c_int64_t) :: t0,t1,rate
        integer :: i,r
        checksum=0
        call system_clock(t0,rate)
        do r=1,repeats
            do i=1,n
                checksum=checksum+rebuilt_q(method(i),a(i),b(i),nu(i),int(setting(i),c_int64_t))
            end do
        end do
        call system_clock(t1)
        sink=checksum
        seconds=real(t1-t0,c_double)/real(rate,c_double)/repeats
    end function

    subroutine screen(n,repeats,a,b,nu,old,new,old_n,new_n,ratio) bind(C)
        use rebuild_adapter, only: rebuilt_time
        integer(c_int), value :: n,repeats
        real(c_double), intent(in) :: a(n),b(n),nu(n)
        integer(c_int), intent(in) :: old(n),new(n),old_n(n),new_n(n)
        real(c_double), intent(out) :: ratio(n)
        integer :: i,r
        real(c_double) :: t_old,t_new,ratios(3)
        integer(c_int64_t) :: count,on,nn
        count=int(repeats,c_int64_t)
        do i=1,n
            on=int(old_n(i),c_int64_t)
            nn=int(new_n(i),c_int64_t)
            do r=1,3
                if (mod(i+r,2)==0) then
                    t_old=rebuilt_time(old(i),a(i),b(i),nu(i),on,count)
                    t_new=rebuilt_time(new(i),a(i),b(i),nu(i),nn,count)
                else
                    t_new=rebuilt_time(new(i),a(i),b(i),nu(i),nn,count)
                    t_old=rebuilt_time(old(i),a(i),b(i),nu(i),on,count)
                end if
                ratios(r)=0
                if (t_new>0) ratios(r)=t_old/t_new
            end do
            ratio(i)=sum(ratios)-minval(ratios)-maxval(ratios)
        end do
    end subroutine
end module
