subroutine xerbla(srname, info)

  implicit none

!  -- LAPACK auxiliary routine (preliminary version) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..

!     .. Scalar Arguments ..
  character(len=*), intent(in) :: srname
  integer, intent(in) :: info

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

  write(*,fmt=9999) srname( 1:len_trim( srname ) ), info
  stop

9999 format(' ** On entry to ', A, ' parameter number ', I2, ' had an illegal value' )

end subroutine xerbla
