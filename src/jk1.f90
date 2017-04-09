!  Pointwise J-Kogbetliantz
module jk1
  implicit none

  include 'params.f90'
  include 'c_ifcs.f90'

contains

  include 'dkogul.f90'
  include 'dkogt2.f90'
  include 'dkogh2.f90'
  include 'jkstep.f90'

  subroutine djksvd(uplo, m, n, G, ldg, U, ldu, V, ldv, jvec, nplus, tol, nsweep, iflags, info)
    implicit none

    !   nsweep .eq. 0: iterate until convergence

    character(len=1), intent(in) :: uplo
    integer, intent(in) :: m, n, ldg, ldu, ldv, jvec(n), nplus, nsweep, iflags
    double precision, intent(in) :: tol
    double precision, intent(inout) :: G(ldg, n), U(ldu, m), V(ldv, n)
    integer, intent(out) :: info(2)

    integer :: maxcyc, sweep, swhalf, step
    integer :: i, j, status
    logical :: upper, rowcyc, hyp

    external :: dscal, dswap
    logical, external :: lsame

    if (lsame(uplo, 'U')) then
       upper = .true.
    else if (lsame(uplo, 'L')) then
       upper = .false.
    else
       info(1) = -1
       info(2) = ichar(uplo)
       return
    end if

    if (m .lt. 0) then
       info(1) = -2
       info(2) = m
       return
    end if

    if ((n .lt. 0) .or. (n .gt. m)) then
       info(1) = -3
       info(2) = n
       return
    end if
    
    if (ldg .lt. m) then
       info(1) = -5
       info(2) = ldg
       return
    end if
    
    if (ldu .lt. m) then
       info(1) = -7
       info(2) = ldu
       return
    end if

    if (ldv .lt. n) then
       info(1) = -9
       info(2) = ldv
       return
    end if

    if ((nplus .lt. 0) .or. (nplus .gt. n)) then
       info(1) = -11
       info(2) = nplus
       return
    end if

    if (tol .lt. ZERO) then
       info(1) = -12
       info(2) = -1
       return
    end if

    if (nsweep .lt. 0) then
       info(1) = -13
       info(2) = nsweep
       return
    end if

    if ((iflags .lt. 0) .or. (iflags .ge. IFLAG_MAXFLG)) then
       info(1) = -14
       info(2) = iflags
       return
    end if

    info = 0

    if (m .eq. 0) return

    rowcyc = (iand(iflags, IFLAG_ROWCYC) .ne. 0)
    if (nsweep .eq. 0) then
       maxcyc = huge(maxcyc)
    else
       maxcyc = nsweep
    end if

    if (rowcyc) then
       do sweep = 1, maxcyc
          do swhalf = 0, 1
             step = 0
             status = 0

             write (*,9) 'Sweep: ', sweep, ' swhalf: ', swhalf

             do i = 1, n - 1
                do j = i + 1, n
                   call jkstep(upper, i, j, m, n, G, ldg, U, ldu, V, ldv, jvec, nplus, tol, sweep, swhalf, step, hyp, info)
                   if (info(1) .ne. 0) return
                   status = max(status, info(2))
                end do
             end do
             if (status .eq. 0) goto 1
             upper = .not. upper
          end do
       end do
    else
       do sweep = 1, maxcyc
          do swhalf = 0, 1
             step = 0
             status = 0
             
             write (*,9) 'Sweep: ', sweep, ' swhalf: ', swhalf

             do j = 2, n
                do i = 1, j - 1
                   call jkstep(upper, i, j, m, n, G, ldg, U, ldu, V, ldv, jvec, nplus, tol, sweep, swhalf, step, hyp, info)
                   if (info(1) .ne. 0) return
                   status = max(status, info(2))
                end do
             end do
             if (status .eq. 0) goto 1
             upper = .not. upper
          end do
       end do
    end if

1   if (sweep .le. maxcyc) then
       info(2) = 2 * sweep + swhalf - 1
    else
       info(2) = -maxcyc
    end if

    if (iand(iflags, IFLAG_PPROCU) .ne. 0) then
       do j = 1, m - 1
          call dswap(m - j, U(j + 1, j), 1, U(j, j + 1), ldu)
       end do
    end if

    if (iand(iflags, IFLAG_PPROCV) .ne. 0) then
       do j = 1, n
          if (jvec(j) .ne. 1) call dscal(n, dble(jvec(j)), V(1, j), 1)
       end do
       do i = 1, n
          if (jvec(i) .ne. 1) call dscal(n, dble(jvec(i)), V(i, 1), ldv)
       end do
       do j = 1, n - 1
          call dswap(n - j, V(j + 1, j), 1, V(j, j + 1), ldv)
       end do
    end if

9   format (A,I2.2,A,I1.1)
  end subroutine djksvd
end module jk1
