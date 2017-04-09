program jk
#ifndef NDEBUG
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: ieee_features
#endif
  use h5data
  use jk1

  implicit none

  integer, parameter :: arglen = 256

  character(len=arglen) :: ah5f, ah5g, ah5r
  integer :: iflags

  integer(hid_t) :: fid, gid

  integer :: idadim(idalen)
  integer :: ldg, n, nrank, nplus
  equivalence (idadim(1), ldg), (idadim(2), n), (idadim(3), nrank), (idadim(4), nplus)
  integer :: ldu, ldv

  double precision, dimension(:,:), allocatable :: G, Q, U, V
  integer, dimension(:), allocatable :: jvec, ipl, invp
  double precision, dimension(:), allocatable :: tau, work

  double precision :: work1(1)
  integer :: lwork, lenwrk, i

  character(len=1) :: uplo
  double precision :: tol
  integer :: nsweep, status

  double precision :: Tstart, Tstop, Tcpu(timlen)
  logical :: fexist
  integer :: info(2)

  external :: dgeqrf, dgeqlf
  external :: dorgqr, dorgql
  external :: dlacpy, dlaset
  external :: dgeqrfp, dscal

  ! EXECUTABLE STATEMENTS

#ifndef NDEBUG
  call ieee_set_halting_mode(ieee_usual, .true.)
#endif

  Tcpu = ZERO

  call readcl(ah5f, ah5g, ah5r, iflags, info)
  if (info(1) .ne. 0) then
     print *, info(1), info(2)
     stop 'Error reading the command line!'
  end if

  inquire(file=ah5f, exist=fexist)
  if (.not. fexist) stop 'Input file does not exist!'

  call h5open_f(info(1))
  if (info(1) .ne. 0) stop 'HDF5 not initialized!'

  call h5fopen_f(ah5f, H5F_ACC_RDONLY_F, fid, info(1))
  if (info(1) .ne. 0) stop 'Error opening the input file!'

  call h5gopen_f(fid, ah5g, gid, info(1))
  if (info(1) .ne. 0) stop 'Error opening the input group!'

  call readh0(gid, idadim, info)
  if (info(1) .ne. 0) stop 'Error reading IDADIM!'

  if ((nplus .lt. 0) .or. (nrank .le. 0) .or. (n .le. 0) .or. (ldg .le. 0)) stop 'IDADIM invalid!'
  if ((nplus .gt. nrank) .or. (nrank .gt. n) .or. (n .gt. ldg)) stop 'IDADIM inconsistent!'

  call dmkldm(n, ldu) ! ldu >= n
  call dmkldm(nrank, ldv) ! ldv >= nrank

  allocate(G(ldg, nrank))
  allocate(jvec(nrank))

  call readh1(gid, n, nrank, G, ldg, jvec, info)
  if (info(1) .ne. 0) then
     print *, info(2)
     stop 'Error reading the input data!'
  end if
  if (ldg .gt. n) call dlaset('A', ldg - n, nrank, ZERO, ZERO, G(n + 1, 1), ldg)

  if (iand(iflags, IFLAG_INVGJP) .ne. 0) then
     allocate(ipl(nrank))
     allocate(invp(nrank))

     call readh2(gid, nrank, ipl, invp, info)
     if (info(1) .ne. 0) then
        print *, info(2)
        stop 'Error reading the permutation data!'
     end if
  end if

  call h5gclose_f(gid, info(1))
  if (info(1) .ne. 0) stop 'Error closing the input group!'

  call h5fclose_f(fid, info(1))
  if (info(1) .ne. 0) stop 'Error closing the input file!'

  if (iand(iflags, IFLAG_INVGJP) .ne. 0) then
     call cpu_time(Tstart)
     call invgjp(n, nrank, G, ldg, jvec, ipl, invp)
     call cpu_time(Tstop)
     
     deallocate(invp)
     deallocate(ipl)

     Tcpu(1) = Tstop - Tstart
     print *, 'INVGJP=', Tcpu(1), 's.'
  else
     Tcpu(1) = ZERO
  end if

  allocate(Q(ldu, n))
  if (ldu .gt. n) call dlaset('A', ldu - n, n, ZERO, ZERO, Q(n + 1, 1), ldu)

  allocate(tau(max(1, min(n, nrank))))

  lwork = -1
  lenwrk = 1

  if (iand(iflags, IFLAG_DGEQRF) .ne. 0) then
     call dgeqrfp(n, nrank, G, ldg, tau, work1, lwork, info(1))
  else
     call dgeqlf(n, nrank, G, ldg, tau, work1, lwork, info(1))
  end if
  if (info(1) .ne. 0) then
     print *, info(1)
     if (iand(iflags, IFLAG_DGEQRF) .ne. 0) then
        stop 'Error in workspace query for DGEQRFP!'
     else
        stop 'Error in workspace query for DGEQLF!'
     end if
  end if

  lenwrk = max(lenwrk, max(1, ceiling(work1(1))))

  if (iand(iflags, IFLAG_DGEQRF) .ne. 0) then
     call dorgqr(n, n, nrank, Q, ldu, tau, work1, lwork, info(1))
  else
     call dorgql(n, n, nrank, Q, ldu, tau, work1, lwork, info(1))
  end if
  if (info(1) .ne. 0) then
     print *, info(1)
     if (iand(iflags, IFLAG_DGEQRF) .ne. 0) then
        stop 'Error in workspace query for DORGQR!'
     else
        stop 'Error in workspace query for DORGQL!'
     end if
  end if

  lenwrk = max(lenwrk, max(1, ceiling(work1(1))))

  lwork = lenwrk
  allocate(work(lwork))

  call cpu_time(Tstart)
  if (iand(iflags, IFLAG_DGEQRF) .ne. 0) then
     call dgeqrfp(n, nrank, G, ldg, tau, work, lwork, info(1))
  else
     call dgeqlf(n, nrank, G, ldg, tau, work, lwork, info(1))
  end if
  if (info(1) .ne. 0) then
     print *, info(1)
     if (iand(iflags, IFLAG_DGEQRF) .ne. 0) then
        stop 'Error in DGEQRFP!'
     else
        stop 'Error in DGEQLF!'
     end if
  end if
  call cpu_time(Tstop)
  Tcpu(2) = Tstop - Tstart
  if (iand(iflags, IFLAG_DGEQRF) .ne. 0) then
     print *, 'DGEQRFP=', Tcpu(2), 's.'
  else
     print *, 'DGEQLF=', Tcpu(2), 's.'
  end if

  call dlacpy('A', n, nrank, G, ldg, Q, ldu)

  call cpu_time(Tstart)
  if (iand(iflags, IFLAG_DGEQRF) .ne. 0) then
     call dorgqr(n, n, nrank, Q, ldu, tau, work, lwork, info(1))
  else
     call dorgql(n, n, nrank, Q, ldu, tau, work, lwork, info(1))
  end if
  if (info(1) .ne. 0) then
     print *, info(1)
     if (iand(iflags, IFLAG_DGEQRF) .ne. 0) then
        stop 'Error in DORGQR!'
     else
        stop 'Error in DORGQL!'
     end if
  end if
  call cpu_time(Tstop)
  Tcpu(3) = Tstop - Tstart
  if (iand(iflags, IFLAG_DGEQRF) .ne. 0) then
     print *, 'DORGQR=', Tcpu(3), 's.'
  else
     print *, 'DORGQL=', Tcpu(3), 's.'
  end if

  deallocate(work)
  deallocate(tau)

  if (iand(iflags, IFLAG_DGEQRF) .ne. 0) then
     uplo = 'U'
     call dlaset('L', n - 1, nrank - 1, ZERO, ZERO, G(2, 1), ldg)
  else
     uplo = 'L'
     call dlaset('U', n - 1, nrank - 1, ZERO, ZERO, G(1, 2), ldg)
  end if

  do i = 1, nrank
     if (G(i, i) .eq. ZERO) then
        stop 'I refuse to work on a singular matrix!'
     else if (G(i, i) .lt. ZERO) then
        call dscal(n, MONE, Q(1, i), 1)
        if (iand(iflags, IFLAG_DGEQRF) .ne. 0) then
           call dscal(nrank - i + 1, MONE, G(i, i), ldg)
        else
           call dscal(i, MONE, G(i, 1), ldg)
        end if
     end if
  end do

  allocate(U(ldu, n))
  if (ldu .gt. n) call dlaset('A', ldu - n, n, ZERO, ZERO, U(n + 1, 1), ldu)

  call dlaset('A', ldu, n, ZERO, ONE, U, ldu) ! U = I

  allocate(V(ldv, nrank))
  if (ldv .gt. nrank) call dlaset('A', ldv - nrank, nrank, ZERO, ZERO, V(nrank + 1, 1), ldv)

  call dlaset('A', ldv, nrank, ZERO, ONE, V, ldv) ! V = I

  inquire(file=ah5r, exist=fexist)
  if (fexist) then
     call h5fopen_f(ah5r, H5F_ACC_RDWR_F, fid, info(1))
     if (info(1) .ne. 0) stop 'Error opening the output file!'
  else
     call h5fcreate_f(ah5r, H5F_ACC_TRUNC_F, fid, info(1))
     if (info(1) .ne. 0) stop 'Error creating the output file!'
  end if

  tol = epsilon(ONE)
  nsweep = 30
  call cpu_time(Tstart)
  call djksvd(uplo, n, nrank, G, ldg, U, ldu, V, ldv, jvec, nplus, tol, nsweep, iflags, info)
  call cpu_time(Tstop)
  if (info(1) .ne. 0) then
     print *, info(1), info(2)
     print *, 'Error in DJKSVD!'
  end if
  Tcpu(4) = Tstop - Tstart
  print *, 'DJKSVD=', Tcpu(4), 's.'
  status = info(2)
  print *, 'DJKSVD:', status

  call h5gcreate_f(fid, ah5g, gid, info(1))
  if (info(1) .ne. 0) stop 'Error creating the output group!'

  call dumph5(gid, n, nrank, G, ldg, jvec, nplus, U, ldu, V, ldv, Q, ldu, tol, iflags, status, Tcpu, info)
  if (info(1) .ne. 0) then
     print *, info(2)
     stop 'Error dumping output file!'
  end if

  deallocate(V)
  deallocate(U)
  deallocate(Q)
  deallocate(jvec)
  deallocate(G)

  call h5gclose_f(gid, info(1))
  if (info(1) .ne. 0) stop 'Error closing the output group!'

  call h5fclose_f(fid, info(1))
  if (info(1) .ne. 0) stop 'Error closing the output file!'

  call h5close_f(info(1))
  if (info(1) .ne. 0) stop 'Error shutting down HDF5!'

#ifndef NDEBUG
  stop 'jk.exe successfully terminated.'
#endif

contains

  include 'dmkldm.f90'
  include 'invgjp.f90'

  subroutine readcl(ah5f, ah5g, ah5r, iflags, info)
    implicit none

    integer, parameter :: argcnt = 4

    character(len=arglen), intent(out) :: ah5f, ah5g, ah5r
    integer, intent(out) :: iflags, info(2)

    character(len=arglen) :: argstr
    integer :: length

    info(1) = command_argument_count()
    if ((info(1) .lt. 0) .or. (info(1) .gt. argcnt)) stop 'jk.exe H5F H5G H5R IFLAGS'

    if (info(1) .ge. 1) then
       call get_command_argument(1, argstr, length, info(2))
       if (info(2) .ne. 0) then
          info(1) = -1
          return
       end if
    else
       print *, 'H5F (input file):'
       read *, argstr
    end if
    ah5f = trim(argstr)

    if (info(1) .ge. 2) then
       call get_command_argument(2, argstr, length, info(2))
       if (info(2) .ne. 0) then
          info(1) = -2
          return
       end if
    else
       print *, 'H5G (input/output group):'
       read *, argstr
    end if
    ah5g = trim(argstr)

    if (info(1) .ge. 3) then
       call get_command_argument(3, argstr, length, info(2))
       if (info(2) .ne. 0) then
          info(1) = -3
          return
       end if
    else
       print *, 'H5R (output file):'
       read *, argstr
    end if
    ah5r = trim(argstr)

    if (info(1) .ge. 4) then
       call get_command_argument(4, argstr, length, info(2))
       if (info(2) .ne. 0) then
          info(1) = -4
          return
       end if
    else
       print *, 'IFLAGS:'
       read *, argstr
    end if
    read (argstr, *) iflags

    info(2) = info(1) - argcnt
    info(1) = 0
  end subroutine readcl
end program jk
