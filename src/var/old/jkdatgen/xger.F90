SUBROUTINE XGER(M, N, ALPHA, X, INCX, Y, INCY, A, LDA)
  IMPLICIT NONE
#include "wp.F"

!     .. Scalar Arguments ..
  REAL(KIND=WP), INTENT(IN) :: ALPHA
  INTEGER, INTENT(IN) :: INCX, INCY, LDA, M, N

!     .. Array Arguments ..
  REAL(KIND=WP), INTENT(INOUT) :: A(LDA,*)
  REAL(KIND=WP), INTENT(IN) :: X(*), Y(*)

!  Purpose
!  =======
!
!  XGER   performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Arguments
!  ==========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL(KIND=WP).
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - REAL(KIND=WP) array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - REAL(KIND=WP) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - REAL(KIND=WP) array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  Further Details
!  ===============
!
!  Level 2 Blas routine.
!
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!  =====================================================================

!     .. Parameters ..
  REAL(KIND=WP), PARAMETER :: ZERO = +0.0E+0_WP

!     .. Local Scalars ..
  REAL(KIND=WP) :: TEMP
  INTEGER :: I, INFO, IX, J, JY, KX

!     .. External Subroutines ..
  EXTERNAL :: XERBLA

!     Test the input parameters.

  INFO = 0
  IF (M .LT. 0) THEN
     INFO = 1
  ELSE IF (N .LT. 0) THEN
     INFO = 2
  ELSE IF (INCX .EQ. 0) THEN
     INFO = 5
  ELSE IF (INCY .EQ. 0) THEN
     INFO = 7
  ELSE IF (LDA .LT. MAX(1,M)) THEN
     INFO = 9
  END IF
  IF (INFO .NE. 0) THEN
     CALL XERBLA('XGER  ',INFO)
     RETURN
  END IF

!     Quick return if possible.

  IF ((M .EQ. 0) .OR. (N .EQ. 0) .OR. (ALPHA .EQ. ZERO)) RETURN

!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.

  IF (INCY .GT. 0) THEN
     JY = 1
  ELSE
     JY = 1 - (N-1)*INCY
  END IF
  IF (INCX .EQ. 1) THEN
     DO J = 1, N
        IF (Y(JY) .NE. ZERO) THEN
           TEMP = ALPHA*Y(JY)
           DO I = 1, M
              A(I,J) = A(I,J) + X(I)*TEMP
           END DO
        END IF
        JY = JY + INCY
     END DO
  ELSE
     IF (INCX .GT. 0) THEN
        KX = 1
     ELSE
        KX = 1 - (M-1)*INCX
     END IF
     DO J = 1, N
        IF (Y(JY) .NE. ZERO) THEN
           TEMP = ALPHA*Y(JY)
           IX = KX
           DO I = 1, M
              A(I,J) = A(I,J) + X(IX)*TEMP
              IX = IX + INCX
           END DO
        END IF
        JY = JY + INCY
     END DO
  END IF

END SUBROUTINE XGER
