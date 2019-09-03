SUBROUTINE XLACPY(UPLO, M, N, A, LDA, B, LDB)
  IMPLICIT NONE
#include "wp.F"

!  -- LAPACK auxiliary routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

!     .. Scalar Arguments ..
  CHARACTER, INTENT(IN) :: UPLO
  INTEGER, INTENT(IN) :: LDA, LDB, M, N

!     .. Array Arguments ..
  REAL(KIND=WP), INTENT(IN) :: A(LDA,*)
  REAL(KIND=WP), INTENT(OUT) :: B(LDB,*)

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
!  A       (input) REAL(KIND=WP) array, dimension (LDA,N)
!          The m by n matrix A.  If UPLO = 'U', only the upper triangle
!          or trapezoid is accessed; if UPLO = 'L', only the lower
!          triangle or trapezoid is accessed.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  B       (output) REAL(KIND=WP) array, dimension (LDB,N)
!          On exit, B = A in the locations specified by UPLO.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,M).
!
!  =====================================================================

!     .. Local Scalars ..
  INTEGER :: I, J

!     .. External Functions ..
  LOGICAL, EXTERNAL :: LSAME

!    .. Executable Statements ..

  IF (LSAME(UPLO,'U')) THEN
     DO J = 1, N
        DO I = 1, MIN(J,M)
           B(I,J) = A(I,J)
        END DO
     END DO
  ELSE IF (LSAME(UPLO,'L')) THEN
     DO J = 1, N
        DO I = J, M
           B(I,J) = A(I,J)
        END DO
     END DO
  ELSE
     DO J = 1, N
        DO I = 1, M
           B(I,J) = A(I,J)
        END DO
     END DO
  END IF

END SUBROUTINE XLACPY
