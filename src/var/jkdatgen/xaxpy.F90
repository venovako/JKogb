SUBROUTINE XAXPY(N, XA, XX, INCX, XY, INCY)
  IMPLICIT NONE
#include "wp.F"

!     .. Scalar Arguments ..
  REAL(KIND=WP), INTENT(IN) :: XA
  INTEGER, INTENT(IN) :: INCX, INCY, N

!     .. Array Arguments ..
  REAL(KIND=WP), INTENT(IN) :: XX(*)
  REAL(KIND=WP), INTENT(INOUT) :: XY(*)

!  Purpose
!  =======
!
!     XAXPY constant times a vector plus a vector.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack.
!
!  =====================================================================

  REAL(KIND=WP), PARAMETER :: ZERO = +0.0E+0_WP

!     .. Local Scalars ..
  INTEGER :: I, IX, IY

  IF (N .LE. 0) RETURN
  IF (XA .EQ. ZERO) RETURN

  IX = 1
  IY = 1
  IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
  IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
  DO I = 1, N
     XY(IY) = XY(IY) + XA*XX(IX)
     IX = IX + INCX
     IY = IY + INCY
  END DO

END SUBROUTINE XAXPY
