SUBROUTINE XSYR(UPLO, N, ALPHA, X, INCX, A, LDA)
  IMPLICIT NONE
#include "wp.F"

!     .. Scalar Arguments ..
  REAL(KIND=WP), INTENT(IN) :: ALPHA
  INTEGER, INTENT(IN) :: INCX, LDA, N
  CHARACTER, INTENT(IN) :: UPLO

!     .. Array Arguments ..
  REAL(KIND=WP), INTENT(INOUT) :: A(LDA,*)
  REAL(KIND=WP), INTENT(IN) :: X(*)

!  Purpose
!  =======
!
!  XSYR   performs the symmetric rank 1 operation
!
!     A := alpha*x*x' + A,
!
!  where alpha is a real scalar, x is an n element vector and A is an
!  n by n symmetric matrix.
!
!  Arguments
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the array A is to be referenced as
!           follows:
!
!              UPLO = 'U' or 'u'   Only the upper triangular part of A
!                                  is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the lower triangular part of A
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL(KIND=WP).
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - REAL(KIND=WP) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  A      - REAL(KIND=WP) array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the symmetric matrix and the strictly
!           lower triangular part of A is not referenced. On exit, the
!           upper triangular part of the array A is overwritten by the
!           upper triangular part of the updated matrix.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the symmetric matrix and the strictly
!           upper triangular part of A is not referenced. On exit, the
!           lower triangular part of the array A is overwritten by the
!           lower triangular part of the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
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
  INTEGER :: I, INFO, IX, J, JX, KX

!     .. External Functions ..
  LOGICAL, EXTERNAL :: LSAME

!     .. External Subroutines ..
  EXTERNAL :: XERBLA

!     Test the input parameters.

  INFO = 0
  IF (.NOT. LSAME(UPLO,'U') .AND. .NOT. LSAME(UPLO,'L')) THEN
     INFO = 1
  ELSE IF (N .LT. 0) THEN
     INFO = 2
  ELSE IF (INCX .EQ. 0) THEN
     INFO = 5
  ELSE IF (LDA .LT. MAX(1,N)) THEN
     INFO = 7
  END IF
  IF (INFO .NE. 0) THEN
     CALL XERBLA('XSYR  ',INFO)
     RETURN
  END IF

!     Quick return if possible.

  IF ((N .EQ. 0) .OR. (ALPHA .EQ. ZERO)) RETURN

!     Set the start point in X if the increment is not unity.

  IF (INCX .LE. 0) THEN
     KX = 1 - (N-1)*INCX
  ELSE IF (INCX .NE. 1) THEN
     KX = 1
  END IF

!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.

  IF (LSAME(UPLO,'U')) THEN

!        Form  A  when A is stored in upper triangle.

     IF (INCX .EQ. 1) THEN
        DO J = 1, N
           IF (X(J) .NE. ZERO) THEN
              TEMP = ALPHA * X(J)
              DO I = 1, J
                 A(I,J) = A(I,J) + X(I) * TEMP
              END DO
           END IF
        END DO
     ELSE
        JX = KX
        DO J = 1, N
           IF (X(JX) .NE. ZERO) THEN
              TEMP = ALPHA * X(JX)
              IX = KX
              DO I = 1, J
                 A(I,J) = A(I,J) + X(IX) * TEMP
                 IX = IX + INCX
              END DO
           END IF
           JX = JX + INCX
        END DO
     END IF
  ELSE

!        Form  A  when A is stored in lower triangle.

     IF (INCX .EQ. 1) THEN
        DO J = 1, N
           IF (X(J) .NE. ZERO) THEN
              TEMP = ALPHA * X(J)
              DO I = J, N
                 A(I,J) = A(I,J) + X(I) * TEMP
              END DO
           END IF
        END DO
     ELSE
        JX = KX
        DO J = 1, N
           IF (X(JX) .NE. ZERO) THEN
              TEMP = ALPHA * X(JX)
              IX = JX
              DO I = J, N
                 A(I,J) = A(I,J) + X(IX) * TEMP
                 IX = IX + INCX
              END DO
           END IF
           JX = JX + INCX
        END DO
     END IF
  END IF

END SUBROUTINE XSYR
