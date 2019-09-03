FUNCTION XNRM2(N, X, INCX)
  IMPLICIT NONE
#include "wp.F"

!     .. Scalar Arguments ..
  INTEGER, INTENT(IN) :: INCX, N

!     .. Array Arguments ..
  REAL(KIND=WP), INTENT(IN) :: X(*)

  REAL(KIND=WP) :: XNRM2

!  Purpose
!  =======
!
!  XNRM2 returns the euclidean norm of a vector via the function
!  name, so that
!
!     XNRM2 := sqrt( x'*x )
!
!  Further Details
!  ===============
!
!     Sven Hammarling, Nag Ltd.
!
!  =====================================================================

!     .. Parameters ..
  REAL(KIND=WP), PARAMETER :: ONE = +1.0E+0_WP, ZERO = +0.0E+0_WP

!     .. Local Scalars ..
  REAL(KIND=WP) :: ABSXI, NORM, SCAL, SSQ
  INTEGER :: IX

  IF (N .LT. 1 .OR. INCX .LT. 1) THEN
     NORM = ZERO
  ELSE IF (N .EQ. 1) THEN
     NORM = ABS(X(1))
  ELSE
     SCAL = ZERO
     SSQ = ONE

!        The following loop is equivalent to this call to the LAPACK
!        auxiliary routine:
!        CALL XLASSQ( N, X, INCX, SCALE, SSQ )

     DO IX = 1, 1 + (N-1)*INCX, INCX
        IF (X(IX) .NE. ZERO) THEN
           ABSXI = ABS(X(IX))
           IF (SCAL .LT. ABSXI) THEN
              SSQ = ONE + SSQ * (SCAL/ABSXI)**2
              SCAL = ABSXI
           ELSE
              SSQ = SSQ + (ABSXI/SCAL)**2
           END IF
        END IF
     END DO
     NORM = SCAL * SQRT(SSQ)
  END IF

  XNRM2 = NORM

END FUNCTION XNRM2
