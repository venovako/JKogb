FUNCTION XDOT(N, XX, INCX, XY, INCY)
  IMPLICIT NONE
#include "wp.F"

!     .. Scalar Arguments ..
  INTEGER, INTENT(IN) :: INCX, INCY, N

!     .. Array Arguments ..
  REAL(KIND=WP), INTENT(IN) :: XX(*), XY(*)

  REAL(KIND=WP) :: XDOT

!  Purpose
!  =======
!
!     XDOT forms the dot product of two vectors.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack.
!
!  =====================================================================

  REAL(KIND=WP), PARAMETER :: ZERO = +0.0E+0_WP

!     .. Local Scalars ..
  REAL(KIND=WP) :: XTEMP
  INTEGER :: I, IX, IY

  XDOT = ZERO
  XTEMP = ZERO

  IF (N .LE. 0) RETURN

  IX = 1
  IY = 1
  IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
  IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
  DO I = 1, N
     XTEMP = XTEMP + XX(IX)*XY(IY)
     IX = IX + INCX
     IY = IY + INCY
  END DO
  XDOT = XTEMP

END FUNCTION XDOT
