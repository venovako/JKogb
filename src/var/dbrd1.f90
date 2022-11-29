!>@brief \b DBRD1 reads a one-dimensional double precision real array from a binary file.
!!
!!@param U [IN]; a connected unit.
!!@param M [IN]; the length of \f$G\f$.
!!@param G [OUT]; a double precision real array to be read.
!!@param INFO [OUT]; zero on success, \f$-i\f$ if the \f$i\f$th argument had an illegal value, or a positive I/O error code.
SUBROUTINE DBRD1(U, M, G, INFO)
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: REAL64
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: U, M
  REAL(KIND=REAL64), INTENT(OUT) :: G(M)
  INTEGER, INTENT(OUT) :: INFO
  INFO = 0
  IF (M .LT. 0) INFO = -2
  IF (INFO .NE. 0) RETURN
  IF (M .EQ. 0) RETURN
  READ (UNIT=U, IOSTAT=INFO) G
END SUBROUTINE DBRD1
