SUBROUTINE XSWAP(N, XX, INCX, XY, INCY)
  IMPLICIT NONE
#include "wp.F"

!     .. Scalar Arguments ..
  INTEGER, INTENT(IN) :: INCX, INCY, N

!     .. Array Arguments ..
  REAL(KIND=WP), INTENT(INOUT) :: XX(*), XY(*)

!  Purpose
!  =======
!
!     interchanges two vectors.
!     uses unrolled loops for increments equal one.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack.
!
!  =====================================================================

!     .. Local Scalars ..
  REAL(KIND=WP) :: XTEMP
  INTEGER :: I, IX, IY

  IF (N .LE. 0) RETURN

  IX = 1
  IY = 1
  IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
  IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
  DO I = 1, N
     XTEMP = XX(IX)
     XX(IX) = XY(IY)
     XY(IY) = XTEMP
     IX = IX + INCX
     IY = IY + INCY
  END DO

END SUBROUTINE XSWAP
