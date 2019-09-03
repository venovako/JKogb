SUBROUTINE XSCAL(N, XA, XX, INCX)
  IMPLICIT NONE
#include "wp.F"

!     .. Scalar Arguments ..
  REAL(KIND=WP), INTENT(IN) :: XA
  INTEGER, INTENT(IN) :: INCX, N

!     .. Array Arguments ..
  REAL(KIND=WP), INTENT(INOUT) :: XX(*)

!  Purpose
!  =======
!
!     XSCAL scales a vector by a constant.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack.
!
!  =====================================================================

!     .. Local Scalars ..
  INTEGER :: I, NINCX

  IF (N .LE. 0 .OR. INCX .LE. 0) RETURN

  NINCX = N*INCX
  DO I = 1, NINCX, INCX
     XX(I) = XA*XX(I)
  END DO

END SUBROUTINE XSCAL
