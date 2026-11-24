! Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
! Paper: Efficient numerical computation of the generalized Marcum
! Q-function with applications (2026).
! Preprint: https://arxiv.org/abs/2611.29209
! Citation: CITATION.cff / CITATION.bib at the repository root.
! Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

program validate
  use iso_c_binding, only: dp => c_double, c_int64_t
  use iso_fortran_env, only: output_unit
  use ieee_arithmetic, only: ieee_is_finite
  use omp_lib, only: omp_get_wtime, omp_get_max_threads
  use MarcumQ, only: marcum
  use IncgamFI, only: incgam
  use paper_algorithm1, only: paper_algorithm1_q
  implicit none

  integer, parameter :: chunk_size = 250000, groups = 4, fields = 10
  integer, parameter :: drawn=1, admitted=2, normal=3, underflow=4, invalid=5, &
                        passed=6, failed=7, route_error=8, passed_underflow=9, flagged_large_bound=10
  integer(c_int64_t), parameter :: max_order = 1000
  integer(c_int64_t) :: count, start, n, totals(fields,groups), batch(fields,groups)
  integer :: length, i, j, seed_size, g, ierr, route_ierr, outcome
  integer :: examples(3,groups)
  integer, allocatable :: seed(:)
  real(dp), allocatable :: draws(:,:)
  real(dp), parameter :: tol=1.0e-10_dp
  real(dp) :: lower(3), span(3), a, b, nu, x, y, log_m, m, pg, qg, p, q, estimate, err, best, &
              max_error(groups), batch_error(groups), t0
  character(32) :: argument

  call get_command_argument(1,argument)
  read(argument,*) count
  if (count < 1) stop 1
  call random_seed(size=seed_size)
  allocate(seed(seed_size),draws(3,chunk_size))
  seed=[(20261109 + 104729*(j-1),j=1,seed_size)]
  call random_seed(put=seed)
  lower=log([0.01_dp,0.001_dp,0.1_dp])
  span=log([100.0_dp,7.0_dp,1000.0_dp])-lower
  totals=0
  max_error=0.0_dp
  examples=0
  write(*,'(a,i0)') 'THREADS ',omp_get_max_threads()
  write(*,'(a,*(i0,1x))') 'SEED ',seed
  write(*,'(a)') 'FIELDS drawn admitted normal underflow invalid passed failed route_error passed_underflow flagged_large_bound'
  write(*,'(a)') 'GROUPS [0.1,1) [1,10) [10,100] (100,1000]'
  t0=omp_get_wtime()

  do start=0_c_int64_t,count-1_c_int64_t,int(chunk_size,c_int64_t)
    length=int(min(int(chunk_size,c_int64_t),count-start))
    call random_number(draws(:,:length))
    batch=0
    batch_error=0.0_dp
    !$omp parallel do schedule(dynamic,256) reduction(+:batch) reduction(max:batch_error) &
    !$omp private(a,b,nu,x,y,g,log_m,m,pg,qg,route_ierr,p,q,ierr,n,estimate,err,best,outcome)
    do i=1,length
      a=exp(lower(1)+span(1)*draws(1,i))
      b=exp(lower(2)+span(2)*draws(2,i))
      nu=exp(lower(3)+span(3)*draws(3,i))
      if (nu < 1.0_dp) then
        g=1
      else if (nu < 10.0_dp) then
        g=2
      else if (nu <= 100.0_dp) then
        g=3
      else
        g=4
      end if
      batch(drawn,g)=batch(drawn,g)+1

      x=a*a/2.0_dp
      y=b*b/2.0_dp
      log_m=log(2.0_dp)-x/2.0_dp+y+nu*log(y)-log_gamma(nu+1.0_dp)
      if (log_m >= log(1.0_dp-2.0e-3_dp)) then
        ! The second bound is at least 2*(1+M) since Q_gamma <= 1.
        if (log_m > log(499.0_dp)) cycle
        m=exp(log_m)
        call incgam(nu,y,pg,qg,route_ierr)
        if (route_ierr /= 0 .or. .not.ieee_is_finite(qg) .or. qg < 0.0_dp .or. qg > 1.0_dp) then
          batch(route_error,g)=batch(route_error,g)+1
          cycle
        end if
        if (2.0_dp*(1.0_dp+m) > 1000.0_dp*qg) cycle
      end if
      batch(admitted,g)=batch(admitted,g)+1

      ! Direct call to the released Fortran routine, also when nu < 1.
      call marcum(nu,x,y,p,q,ierr)
      outcome=0
      if (ierr == 1 .and. p == 0.0_dp .and. q == 1.0_dp) then
        batch(underflow,g)=batch(underflow,g)+1
        if (log_m > log(1.0e-12_dp)) batch(flagged_large_bound,g)=batch(flagged_large_bound,g)+1
        outcome=1
      else if (ierr /= 0 .or. .not.ieee_is_finite(p) .or. .not.ieee_is_finite(q) .or. &
               q <= 0.0_dp .or. q > 1.0_dp .or. p < 0.0_dp .or. p > 1.0_dp .or. &
               abs(p+q-1.0_dp) > 1.0e-13_dp) then
        batch(invalid,g)=batch(invalid,g)+1
        call record_example(2,g,start+i,a,b,nu,log_m,ierr,p,q,-1.0_dp)
        cycle
      else
        batch(normal,g)=batch(normal,g)+1
      end if

      best=huge(1.0_dp)
      do n=1,max_order
        estimate=paper_algorithm1_q(a,b,nu,n)
        if (.not.ieee_is_finite(estimate)) cycle
        err=abs(estimate-q)/q
        best=min(best,err)
        if (err <= tol) exit
      end do
      if (outcome == 1 .and. start == 0) call record_example(1,g,start+i,a,b,nu,log_m,ierr,p,q,best)
      if (best <= tol) then
        batch(passed,g)=batch(passed,g)+1
        if (outcome == 1) batch(passed_underflow,g)=batch(passed_underflow,g)+1
        batch_error(g)=max(batch_error(g),best)
      else
        batch(failed,g)=batch(failed,g)+1
        call record_example(3,g,start+i,a,b,nu,log_m,ierr,p,q,best)
      end if
    end do
    !$omp end parallel do
    totals=totals+batch
    max_error=max(max_error,batch_error)
    if (mod(start+length,1000000_c_int64_t) == 0 .or. start+length == count) then
      write(*,'(a,i0,1x,f12.3,4(1x,i0))') 'PROGRESS ',start+length,omp_get_wtime()-t0, &
          sum(totals(admitted,:)),sum(totals(passed,:)),sum(totals(failed,:)),sum(totals(invalid,:))
      flush(output_unit)
    end if
  end do
  do g=1,groups
    write(*,'(a,i0,1x,*(i0,1x))') 'RESULT ',g,totals(:,g)
    write(*,'(a,i0,1x,es26.17e3)') 'MAX_PASS_ERROR ',g,max_error(g)
  end do
  if (sum(totals(drawn,:)) /= count) stop 2
  if (any(totals(admitted,:) /= totals(normal,:)+totals(underflow,:)+totals(invalid,:))) stop 3
  if (any(totals(admitted,:) /= totals(passed,:)+totals(failed,:)+totals(invalid,:))) stop 4
  if (any(totals(route_error,:) /= 0)) stop 5
  write(*,'(a,f12.3)') 'COMPLETE ',omp_get_wtime()-t0

contains

  subroutine record_example(kind,group,index,a,b,nu,log_m,ierr,p,q,best)
    integer,intent(in) :: kind,group,ierr
    integer(c_int64_t),intent(in) :: index
    real(dp),intent(in) :: a,b,nu,log_m,p,q,best
    !$omp critical(example_output)
    if (examples(kind,group) < 4) then
      write(*,'(a,4(i0,1x),7(es26.17e3,1x))') 'EXAMPLE ',kind,group,index,ierr,a,b,nu,log_m,p,q,best
      examples(kind,group)=examples(kind,group)+1
    end if
    !$omp end critical(example_output)
  end subroutine record_example
end program validate
