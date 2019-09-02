subroutine xlagsy(n, k, d, a, lda, iseed, work, info)

  implicit none

!  -- LAPACK auxiliary test routine (version 3.1)
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..

!     .. Scalar Arguments ..
  integer, intent(out) :: info
  integer, intent(in) :: k, lda, n

!     .. Array Arguments ..
  integer, intent(inout) :: iseed(4)
  real(kind=16), intent(out) :: a(lda,*)
  real(kind=16), intent(in) :: d(*)
  real(kind=16), intent(out) :: work(*)

!  Purpose
!  =======
!
!  XLAGSY generates a real symmetric matrix A, by pre- and post-
!  multiplying a real diagonal matrix D with a random orthogonal matrix:
!  A = U*D*U'. The semi-bandwidth may then be reduced to k by additional
!  orthogonal transformations.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  K       (input) INTEGER
!          The number of nonzero subdiagonals within the band of A.
!          0 <= K <= N-1.
!
!  D       (input) REAL(KIND=16) array, dimension (N)
!          The diagonal elements of the diagonal matrix D.
!
!  A       (output) REAL(KIND=16) array, dimension (LDA,N)
!          The generated n by n symmetric matrix A (the full matrix is
!          stored).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= N.
!
!  ISEED   (input/output) INTEGER array, dimension (4)
!          On entry, the seed of the random number generator; the array
!          elements must be between 0 and 4095, and ISEED(4) must be
!          odd.
!          On exit, the seed is updated.
!
!  WORK    (workspace) REAL(KIND=16) array, dimension (2*N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================

!     .. Parameters ..
  real(kind=16), parameter :: zero = +0.0e+0_16, one = +1.0e+0_16, half = +0.5e+0_16

!     .. Local Scalars ..
  integer :: i, j
  real(kind=16) :: alpha, tau, wa, wb, wn

!     .. External Subroutines ..
  external :: xaxpy, xgemv, xger, xlarnv, xscal, xsymv, xsyr2, xerbla

!     .. External Functions ..
  real(kind=16), external :: xdot, xnrm2

!     .. Executable Statements ..

!     Test the input arguments

  info = 0
  if (n .lt. 0) then
     info = -1
  else if (k .lt. 0 .or. k .gt. n-1) then
     info = -2
  else if (lda .lt. max(1,n)) then
     info = -5
  end if
  if (info .lt. 0) then
     call xerbla('XLAGSY', -info)
     return
  end if

!     initialize lower triangle of A to diagonal matrix

  do j = 1, n
     do i = j + 1, n
        a( i, j ) = zero
     end do
  end do
  do i = 1, n
     a( i, i ) = d( i )
  end do

!     Generate lower triangle of symmetric matrix

  do i = n - 1, 1, -1

!        generate random reflection

     call xlarnv(3, iseed, n-i+1, work)
     wn = xnrm2(n-i+1, work, 1)
     wa = sign(wn, work( 1 ))
     if (wn .eq. zero) then
        tau = zero
     else
        wb = work( 1 ) + wa
        call xscal(n-i, one / wb, work( 2 ), 1)
        work( 1 ) = one
        tau = wb / wa
     end if

!        apply random reflection to A(i:n,i:n) from the left
!        and the right

!        compute  y := tau * A * u

     call xsymv('Lower', n-i+1, tau, a( i, i ), lda, work, 1, zero, work( n+1 ), 1)

!        compute  v := y - 1/2 * tau * ( y, u ) * u

     alpha = -half * tau * xdot( n-i+1, work( n+1 ), 1, work, 1 )
     call xaxpy(n-i+1, alpha, work, 1, work( n+1 ), 1)

!        apply the transformation as a rank-2 update to A(i:n,i:n)

     call xsyr2('Lower', n-i+1, -one, work, 1, work( n+1 ), 1, a( i, i ), lda)
  end do

!     Reduce number of subdiagonals to K

  do i = 1, n - 1 - k

!        generate reflection to annihilate A(k+i+1:n,i)

     wn = xnrm2(n-k-i+1, a( k+i, i ), 1)
     wa = sign(wn, a( k+i, i ))
     if (wn .eq. zero) then
        tau = zero
     else
        wb = a( k+i, i ) + wa
        call xscal(n-k-i, one / wb, a( k+i+1, i ), 1)
        a( k+i, i ) = one
        tau = wb / wa
     end if

!        apply reflection to A(k+i:n,i+1:k+i-1) from the left

     call xgemv('Transpose', n-k-i+1, k-1, one, a( k+i, i+1 ), lda, a( k+i, i ), 1, zero, work, 1)
     call xger(n-k-i+1, k-1, -tau, a( k+i, i ), 1, work, 1, a( k+i, i+1 ), lda)

!        apply reflection to A(k+i:n,k+i:n) from the left and the right

!        compute  y := tau * A * u

     call xsymv('Lower', n-k-i+1, tau, a( k+i, k+i ), lda, a( k+i, i ), 1, zero, work, 1)

!        compute  v := y - 1/2 * tau * ( y, u ) * u

     alpha = -half * tau * xdot(n-k-i+1, work, 1, a( k+i, i ), 1)
     call xaxpy(n-k-i+1, alpha, a( k+i, i ), 1, work, 1)

!        apply symmetric rank-2 update to A(k+i:n,k+i:n)

     call xsyr2('Lower', n-k-i+1, -one, a( k+i, i ), 1, work, 1, a( k+i, k+i ), lda)

     a( k+i, i ) = -wa
     do j = k + i + 1, n
        a( j, i ) = zero
     end do
  end do

!     Store full symmetric matrix

  do j = 1, n
     do i = j + 1, n
        a( j, i ) = a( i, j )
     end do
  end do

end subroutine xlagsy
