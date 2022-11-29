!>@brief \b ZBWR2 writes a two-dimensional double precision complex array to a binary file.
!!
!!@param U [IN]; a connected unit.
!!@param M [IN]; the number of rows of \f$G\f$.
!!@param N [IN]; the number of columns of \f$G\f$.
!!@param G [IN]; a double precision complex array to be written.
!!@param LDG [IN]; the leading dimension of \f$G\f$.
!!@param INFO [OUT]; zero on success, \f$-i\f$ if the \f$i\f$th argument had an illegal value, or a positive I/O error code.
SUBROUTINE ZBWR2(U, M, N, G, LDG, INFO)
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: REAL64
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: U, M, N, LDG
  COMPLEX(KIND=REAL64), INTENT(IN) :: G(LDG,N)
  INTEGER, INTENT(OUT) :: INFO
  INTEGER :: J
  EXTERNAL :: ZBWR1
  INFO = 0
  IF (LDG .LT. MAX(M, 0)) INFO = -5
  IF (N .LT. 0) INFO = -3
  IF (M .LT. 0) INFO = -2
  IF (INFO .NE. 0) RETURN
  IF (M .EQ. 0) RETURN
  IF (N .EQ. 0) RETURN
  DO J = 1, N
     CALL ZBWR1(U, M, G(1,J), INFO)
     IF (INFO .NE. 0) EXIT
  END DO
END SUBROUTINE ZBWR2
