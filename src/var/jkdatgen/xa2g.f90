subroutine xa2g(n, lam, iseed, nrank, a, g, ldg, ipiv, jvec, nplus, ipl, invp, info)

  implicit none

  integer, intent(in) :: n, ldg
  real(kind=8), intent(in) :: lam(n)
  integer, intent(inout) :: iseed(4)
  real(kind=8), intent(out) :: a(ldg,n), g(ldg,n)
  integer, intent(out) :: nrank, ipiv(n), jvec(n), nplus, ipl(n), invp(n), info(2)

  real(kind=8), parameter :: dzero = +0.0e+0_8

  real(kind=16), allocatable :: xlam(:), xa(:,:), work(:)
  integer :: lda, i, j

  external :: jpart, xlacpy, xlagsy, xsybpc

  lda = ldg
  info(1) = 0
  info(2) = 0

  allocate(xa(lda,n))
  allocate(work(2*n))

  allocate(xlam(n))
  do i = 1, n
     xlam(i) = real(lam(i), 16)
  end do
1 call xlagsy(n, n-1, xlam, xa, lda, iseed, work, info(2))
  deallocate(xlam)
  if (info(2) .ne. 0) then
     info(1) = 1
     goto 3
  end if

  do j = 1, n
     do i = 1, n
        a(i,j) = real(xa(i,j), 8)
     end do
  end do
  do j = 1, n
     do i = n + 1, ldg
        a(i,j) = dzero
     end do
  end do

  ipiv = 0
2 call xsybpc('L', n, xa, lda, nrank, work, ipiv, jvec, info(2))
  if (info(2) .eq. 0) then
     call jpart(n, nrank, xa, lda, jvec, nplus, ipl, invp)
     do j = 1, n
        do i = 1, n
           g(i,j) = real(xa(i,j), 8)
        end do
     end do
     do j = 1, n
        do i = n + 1, ldg
           g(i,j) = dzero
        end do
     end do
  else
     info(1) = 2
  end if

3 deallocate(work)
  deallocate(xa)

end subroutine xa2g
