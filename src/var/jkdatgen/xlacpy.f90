subroutine xlacpy(uplo, m, n, a, lda, b, ldb)

  implicit none

!  -- LAPACK auxiliary routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

!     .. Scalar Arguments ..
  character, intent(in) :: uplo
  integer, intent(in) :: lda, ldb, m, n

!     .. Array Arguments ..
  real(kind=16), intent(in) :: a(lda,*)
  real(kind=16), intent(out) :: b(ldb,*)

!  Purpose
!  =======
!
!  XLACPY copies all or part of a two-dimensional matrix A to another
!  matrix B.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies the part of the matrix A to be copied to B.
!          = 'U':      Upper triangular part
!          = 'L':      Lower triangular part
!          Otherwise:  All of the matrix A
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input) REAL(KIND=16) array, dimension (LDA,N)
!          The m by n matrix A.  If UPLO = 'U', only the upper triangle
!          or trapezoid is accessed; if UPLO = 'L', only the lower
!          triangle or trapezoid is accessed.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  B       (output) REAL(KIND=16) array, dimension (LDB,N)
!          On exit, B = A in the locations specified by UPLO.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,M).
!
!  =====================================================================

!     .. Local Scalars ..
  integer :: i, j

!     .. External Functions ..
  logical, external :: lsame

!    .. Executable Statements ..

  if (lsame(uplo,'U')) then
     do j = 1, n
        do i = 1, min(j,m)
           b(i,j) = a(i,j)
        end do
     end do
  else if (lsame(uplo,'L')) then
     do j = 1, n
        do i = j, m
           b(i,j) = a(i,j)
        end do
     end do
  else
     do j = 1, n
        do i = 1, m
           b(i,j) = a(i,j)
        end do
     end do
  end if

end subroutine xlacpy
