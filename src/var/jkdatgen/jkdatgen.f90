program jkdatgen

  use HDF5
  use H5LT

  implicit none

  ! command-line parameters
  integer :: lgix, seed, n, ldg
  character(len=256) :: fil, grp

  logical :: fex
  integer(hid_t) :: fid, gid

  integer :: iseed(4), nrank, nplus, info(2)
  real(kind=8), allocatable :: lam(:), a(:,:), g(:,:)
  integer, allocatable :: ipiv(:), jvec(:), ipl(:), invp(:)

  call readcl(lgix, seed, n, ldg, fil, grp, info)
  if (info(1) .ne. 0) then
    print *, info(1), info(2)
    stop 'readcl'
  end if

  allocate(lam(n))
  allocate(a(ldg,n))
  allocate(g(ldg,n))
  allocate(ipiv(n))
  allocate(jvec(n))
  allocate(ipl(n))
  allocate(invp(n))

  call gendat(lgix, seed, n, lam, iseed, nrank, a, g, ldg, ipiv, jvec, nplus, ipl, invp, info)
  if (info(1) .ne. 0) then
    print *, info(1), info(2)
    stop 'gendat'
  end if

  call h5open_f(info(1))
  if (info(1) .ne. 0) then
    print *, info(1)
    stop 'h5open_f'
  end if

  inquire(file=fil, exist=fex)
  if (fex) then
    call h5fopen_f(fil, H5F_ACC_RDWR_F, fid, info(1))
  else
    call h5fcreate_f(fil, H5F_ACC_TRUNC_F, fid, info(1))
  end if
  if (info(1) .ne. 0) then
    print *, info(1)
    if (fex) then
      stop 'h5open_f'
    else
      stop 'h5create_f'
    end if
  end if

  call h5gcreate_f(fid, grp, gid, info(1))
  if (info(1) .ne. 0) then
    print *, info(1)
    stop 'h5gcreate_f'
  end if

  call h5wrds(gid, n, lam, iseed, a, g, ldg, nrank, ipiv, jvec, nplus, ipl, invp, info)

  deallocate(invp)
  deallocate(ipl)
  deallocate(jvec)
  deallocate(ipiv)
  deallocate(g)
  deallocate(a)
  deallocate(lam)

  if (info(1) .ne. 0) then
    print *, info(1), info(2)
    stop 'h5wrds'
  end if

  call h5gclose_f(gid, info(1))
  if (info(1) .ne. 0) then
    print *, info(1)
    stop 'h5gclose_f'
  end if

  call h5fclose_f(fid, info(1))
  if (info(1) .ne. 0) then
    print *, info(1)
    stop 'h5fclose_f'
  end if

  call h5close_f(info(1))
  if (info(1) .ne. 0) then
    print *, info(1)
    stop 'h5close_f'
  end if

contains

  subroutine readcl(lgix, seed, n, ldg, fil, grp, info)

    implicit none

    integer, intent(out) :: lgix, seed, n, ldg
    character(len=*), intent(out) :: fil, grp
    integer, intent(out) :: info(2)

    character(len=256) :: cas

    info(1) = 0
    info(2) = 0

    if (command_argument_count() .lt. 6) then
      print *, 'jkdatgen.exe LAMGEN SEEDIX N LDG FILE.h5 GROUP [ LAMGEN_PARAMS ]'
      print *, '>> COMMAND LINE (INPUT) ARGUMENTS <<'
      print *, 'LAMGEN : spectrum [\Lambda] generator to use (see below): 1 or 2'
      print *, 'SEEDIX : index of hard-coded pRNG seed (see seedix.f90): 1 or 2'
      print *, 'N      : order of the output matrix: > 0'
      print *, 'LDG    : leading dimension of the output matrix: >= N'
      print *, 'FILE.h5: output HDF5 file (may exist): max 256 chars'
      print *, 'GROUP  : output HDF5 group (must NOT exist): max 256 chars'
      print *, 'LAMGEN | LAMGEN_PARAMS'
      print *, '1      : LAMFILE.txt: max 256 chars, >= N lines [each line = one eigenvalue]'
      print *, '2      : IDIST EPS SCALE'
      print *, ' IDIST : 1 [uniform (0,1)], 2 [uniform(-1,1)], or 3 [normal(\mu=0,\sigma=1)]'
      print *, ' EPS   : pseudorandom \lambda''_i survives iff |\lambda''_i| > EPS'
      print *, ' SCALE : final \lambda_i = \lambda''_i * SCALE'
      print *, '<< OUTPUT DATASETS IN FILE.h5/GROUP >>'
      print *, 'IDADIM : integer(4) { LDG, N, NRANK, NPLUS }'
      print *, ' NRANK : numerical column-rank: <= N [but should fail when rank-deficient]'
      print *, ' NPLUS : number of positive eigenvalues: <= N [+/-0 counts as non-positive]'
      print *, 'LAMBDA : double precision(N): the eigenvalues [prescribed or pseudorandom]'
      print *, 'ISEED  : integer(4): initial seed for (d|x)laran pRNG (see LaPACK dlaran.f)'
      print *, 'A      : double precision(LDG,N): the generated symmetric indefinite matrix'
      print *, 'G      : double precision(LDG,N): the `lower'' Bunch-Parlett factor of'
      print *, '       : P A P^T == G'' Q Q^T J'' Q Q^T G''^T == G J G^T'
      print *, 'IPIV   : integer(N): complete pivoting permutation vector P (see xsybpc.f90)'
      print *, 'JVEC   : integer(N): diag(J) == diag(I_{NPLUS},-I_{NRANK-NPLUS},0_{N-NRANK})'
      print *, 'IPL    : integer(N): Q [J-partitioning permutation vector (see jpart.f90)]'
      print *, 'INVP   : integer(N): Q^T [inverse of IPL (see jpart.f90)]'
      info(1) = -1
      info(2) = command_argument_count() - 6
      return
    end if

    call get_command_argument(1, cas)
    read (cas,*) lgix
    if (lgix .lt. 1) then
      info(1) = 1
      info(2) = lgix
      return
    end if

    call get_command_argument(2, cas)
    read (cas,*) seed
    if (seed .lt. 1) then
      info(1) = 2
      info(2) = seed
      return
    end if

    call get_command_argument(3, cas)
    read (cas,*) n
    if (n .lt. 1) then
      info(1) = 3
      info(2) = n
      return
    end if

    call get_command_argument(4, cas)
    read (cas,*) ldg
    if (ldg .lt. n) then
      info(1) = 4
      info(2) = ldg
      return
    end if

    call get_command_argument(5, fil)
    if (len_trim(fil) .eq. 0) then
      info(1) = 5
      return
    end if

    call get_command_argument(6, grp)
    if (len_trim(grp) .eq. 0) then
      info(1) = 6
      return
    end if

  end subroutine readcl

  subroutine h5wrds(gid, n, lam, iseed, a, g, ldg, nrank, ipiv, jvec, nplus, ipl, invp, info)

    implicit none

    integer(hid_t), intent(in) :: gid
    integer, intent(in) :: n, iseed(4), ldg, nrank, ipiv(n), jvec(n), nplus, ipl(n), invp(n)
    real(kind=8), intent(in) :: lam(n), a(ldg,n), g(ldg,n)
    integer, intent(out) :: info(2)

    integer :: idadim(4)
    integer(hsize_t) :: dims(2)

    info(1) = 0
    info(2) = 0

    idadim = (/ ldg, n, nrank, nplus /)
    dims(1) = 4
1   call h5ltmake_dataset_int_f(gid, 'IDADIM', 1, dims, idadim, info(2))
    if (info(2) .ne. 0) then
      info(1) = 1
      return
    end if

    dims(1) = n
2   call h5ltmake_dataset_double_f(gid, 'LAMBDA', 1, dims, lam, info(2))
    if (info(2) .ne. 0) then
      info(1) = 2
      return
    end if

    dims(1) = 4
3   call h5ltmake_dataset_int_f(gid, 'ISEED', 1, dims, iseed, info(2))
    if (info(2) .ne. 0) then
      info(1) = 3
      return
    end if

    dims(1) = ldg
    dims(2) = n
4   call h5ltmake_dataset_double_f(gid, 'A', 2, dims, a, info(2))
    if (info(2) .ne. 0) then
      info(1) = 4
      return
    end if

    dims(1) = ldg
    dims(2) = n
5   call h5ltmake_dataset_double_f(gid, 'G', 2, dims, g, info(2))
    if (info(2) .ne. 0) then
      info(1) = 5
      return
    end if

    dims(1) = n
6   call h5ltmake_dataset_int_f(gid, 'IPIV', 1, dims, ipiv, info(2))
    if (info(2) .ne. 0) then
      info(1) = 6
      return
    end if

    dims(1) = n
7   call h5ltmake_dataset_int_f(gid, 'JVEC', 1, dims, jvec, info(2))
    if (info(2) .ne. 0) then
      info(1) = 7
      return
    end if

    dims(1) = n
8   call h5ltmake_dataset_int_f(gid, 'IPL', 1, dims, ipl, info(2))
    if (info(2) .ne. 0) then
      info(1) = 8
      return
    end if

    dims(1) = n
9   call h5ltmake_dataset_int_f(gid, 'INVP', 1, dims, invp, info(2))
    if (info(2) .ne. 0) then
      info(1) = 9
      return
    end if

  end subroutine h5wrds

end program jkdatgen
