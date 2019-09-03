SUBROUTINE XERBLA(SRNAME, INFO)
  IMPLICIT NONE

!  -- LAPACK auxiliary routine (preliminary version) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..

!     .. Scalar Arguments ..
  CHARACTER(LEN=*), INTENT(IN) :: SRNAME
  INTEGER, INTENT(IN) :: INFO

!  Purpose
!  =======
!
!  XERBLA  is an error handler for the LAPACK routines.
!  It is called by an LAPACK routine if an input parameter has an
!  invalid value.  A message is printed and execution stops.
!
!  Installers may consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!
!  Arguments
!  =========
!
!  SRNAME  (input) CHARACTER*(*)
!          The name of the routine which called XERBLA.
!
!  INFO    (input) INTEGER
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
!
! =====================================================================

!     .. Executable Statements ..

  WRITE(*,FMT=9999) SRNAME( 1:LEN_TRIM( SRNAME ) ), INFO
  STOP

9999 FORMAT(' ** On entry to ', A, ' parameter number ', I2, ' had an illegal value')

END SUBROUTINE XERBLA
