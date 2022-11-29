!>@brief \b JBRD1 reads a one-dimensional INT64 array from a binary file.
!!
!!@param U [IN]; a connected unit.
!!@param M [IN]; the length of \f$G\f$.
!!@param G [OUT]; an INT64 array to be read.
!!@param INFO [OUT]; zero on success, \f$-i\f$ if the \f$i\f$th argument had an illegal value, or a positive I/O error code.
SUBROUTINE JBRD1(U, M, G, INFO)
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: INT64
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: U, M
  INTEGER(KIND=INT64), INTENT(OUT) :: G(M)
  INTEGER, INTENT(OUT) :: INFO
  INFO = 0
  IF (M .LT. 0) INFO = -2
  IF (INFO .NE. 0) RETURN
  IF (M .EQ. 0) RETURN
  READ (UNIT=U, IOSTAT=INFO) G
END SUBROUTINE JBRD1
