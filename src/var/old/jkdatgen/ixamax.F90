FUNCTION IXAMAX(N, XX, INCX)
  IMPLICIT NONE
#include "wp.F"

!     .. Scalar Arguments ..
  INTEGER, INTENT(IN) :: INCX, N

!     .. Array Arguments ..
  REAL(KIND=WP), INTENT(IN) :: XX(*)

  INTEGER :: IXAMAX

!  Purpose
!  =======
!
!     IXAMAX finds the index of element having max. absolute value.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack.
!
!  =====================================================================

!     .. Local Scalars ..
  REAL(KIND=WP) :: XMAX
  INTEGER :: I, IX

  IXAMAX = 0
  IF (N .LT. 1 .OR. INCX .LE. 0) RETURN
  IXAMAX = 1
  IF (N .EQ. 1) RETURN
  IF (INCX .EQ. 1) GOTO 10

!        code for increment not equal to 1

  IX = 1
  XMAX = ABS(XX(1))
  IX = IX + INCX
  DO I = 2, N
     IF (ABS(XX(IX)) .LE. XMAX) GOTO 5
     IXAMAX = I
     XMAX = ABS(XX(IX))
5    IX = IX + INCX
  END DO
  RETURN

!        code for increment equal to 1

10 XMAX = ABS(XX(1))
  DO I = 2, N
     IF (ABS(XX(I)) .LE. XMAX) CYCLE
     IXAMAX = I
     XMAX = ABS(XX(I))
  END DO

END FUNCTION IXAMAX
