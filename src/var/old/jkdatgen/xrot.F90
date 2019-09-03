SUBROUTINE XROT(N, XX, INCX, XY, INCY, C, S)
  IMPLICIT NONE
#include "wp.F"

!     .. Scalar Arguments ..
  REAL(KIND=WP), INTENT(IN) :: C, S
  INTEGER, INTENT(IN) :: INCX, INCY, N

!     .. Array Arguments ..
  REAL(KIND=WP), INTENT(INOUT) :: XX(*), XY(*)

!  Purpose
!  =======
!
!     XROT applies a plane rotation.
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
  IF (INCX .EQ. 1 .AND. INCY .EQ. 1) GOTO 10

!       code for unequal increments or equal increments not equal
!         to 1

  IX = 1
  IY = 1
  IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
  IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
  DO I = 1, N
     XTEMP = C*XX(IX) + S*XY(IY)
     XY(IY) = C*XY(IY) - S*XX(IX)
     XX(IX) = XTEMP
     IX = IX + INCX
     IY = IY + INCY
  END DO
  RETURN

!       code for both increments equal to 1

10 DO I = 1, N
     XTEMP = C*XX(I) + S*XY(I)
     XY(I) = C*XY(I) - S*XX(I)
     XX(I) = XTEMP
  END DO

END SUBROUTINE XROT
