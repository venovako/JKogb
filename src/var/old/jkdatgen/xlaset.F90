SUBROUTINE XLASET(UPLO, M, N, ALPHA, BETA, A, LDA)
  IMPLICIT NONE
#include "wp.F"

!  -- LAPACK auxiliary routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

!     .. Scalar Arguments ..
  CHARACTER, INTENT(IN) :: UPLO
  INTEGER, INTENT(IN) :: LDA, M, N
  REAL(KIND=WP), INTENT(IN) :: ALPHA, BETA

!     .. Array Arguments ..
  REAL(KIND=WP), INTENT(OUT) :: A(LDA,*)

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
!  ALPHA   (input) REAL(KIND=WP)
!          The constant to which the offdiagonal elements are to be set.
!
!  BETA    (input) REAL(KIND=WP)
!          The constant to which the diagonal elements are to be set.
!
!  A       (input/output) REAL(KIND=WP) array, dimension (LDA,N)
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
  INTEGER :: I, J

!     .. External Functions ..
  LOGICAL, EXTERNAL :: LSAME

!     .. Executable Statements ..

  IF (LSAME(UPLO,'U')) THEN

!        Set the strictly upper triangular or trapezoidal part of the
!        array to ALPHA.

     DO J = 2, N
        DO I = 1, MIN(J-1,M)
           A(I,J) = ALPHA
        END DO
     END DO

  ELSE IF (LSAME(UPLO,'L')) THEN

!        Set the strictly lower triangular or trapezoidal part of the
!        array to ALPHA.

     DO J = 1, MIN(M,N)
        DO I = J+1, M
           A(I,J) = ALPHA
        END DO
     END DO

  ELSE

!        Set the leading m-by-n submatrix to ALPHA.

     DO J = 1, N
        DO I = 1, M
           A(I,J) = ALPHA
        END DO
     END DO
  END IF

!     Set the first min(M,N) diagonal elements to BETA.

  DO I = 1, MIN(M,N)
     A(I,I) = BETA
  END DO

END SUBROUTINE XLASET
