subroutine xlaset(uplo, m, n, alpha, beta, a, lda)

  implicit none

!  -- LAPACK auxiliary routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

!     .. Scalar Arguments ..
  character, intent(in) :: uplo
  integer, intent(in) :: lda, m, n
  real(kind=16), intent(in) :: alpha, beta

!     .. Array Arguments ..
  real(kind=16), intent(out) :: a(lda,*)

!  Purpose
!  =======
!
!  XLASET initializes an m-by-n matrix A to BETA on the diagonal and
!  ALPHA on the offdiagonals.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies the part of the matrix A to be set.
!          = 'U':      Upper triangular part is set; the strictly lower
!                      triangular part of A is not changed.
!          = 'L':      Lower triangular part is set; the strictly upper
!                      triangular part of A is not changed.
!          Otherwise:  All of the matrix A is set.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  ALPHA   (input) REAL(KIND=16)
!          The constant to which the offdiagonal elements are to be set.
!
!  BETA    (input) REAL(KIND=16)
!          The constant to which the diagonal elements are to be set.
!
!  A       (input/output) REAL(KIND=16) array, dimension (LDA,N)
!          On exit, the leading m-by-n submatrix of A is set as follows:
!
!          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
!          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
!          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
!
!          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
! =====================================================================

!     .. Local Scalars ..
  integer :: i, j

!     .. External Functions ..
  logical, external :: lsame

!     .. Executable Statements ..

  if (lsame(uplo,'U')) then

!        Set the strictly upper triangular or trapezoidal part of the
!        array to ALPHA.

     do j = 2, n
        do i = 1, min(j-1,m)
           a(i,j) = alpha
        end do
     end do

  else if (lsame(uplo,'L')) then

!        Set the strictly lower triangular or trapezoidal part of the
!        array to ALPHA.

     do j = 1, min(m,n)
        do i = j+1, m
           a(i,j) = alpha
        end do
     end do

  else

!        Set the leading m-by-n submatrix to ALPHA.

     do j = 1, n
        do i = 1, m
           a(i,j) = alpha
        end do
     end do
  end if

!     Set the first min(M,N) diagonal elements to BETA.

  do i = 1, min(m,n)
     a(i,i) = beta
  end do

end subroutine xlaset
