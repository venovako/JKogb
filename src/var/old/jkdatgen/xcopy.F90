SUBROUTINE XCOPY(N, XX, INCX, XY, INCY)
  IMPLICIT NONE
#include "wp.F"

!     .. Scalar Arguments ..
  INTEGER, INTENT(IN) :: INCX, INCY, N

!     .. Array Arguments ..
  REAL(KIND=WP), INTENT(IN) :: XX(*)
  REAL(KIND=WP), INTENT(OUT) :: XY(*)

!  Purpose
!  =======
!
!     XCOPY copies a vector, x, to a vector, y.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack.
!
!  =====================================================================

!     .. Local Scalars ..
  INTEGER :: I, IX, IY

  IF (N .LE. 0) RETURN

!        code for unequal increments or equal increments
!          not equal to 1

  IX = 1
  IY = 1
  IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
  IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
  DO I = 1, N
     XY(IY) = XX(IX)
     IX = IX + INCX
     IY = IY + INCY
  END DO

END SUBROUTINE XCOPY
