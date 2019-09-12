MODULE DTRANSF
  USE PARAMS
  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE DTRANSFA(A, J, U, Z)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: A(2,2)
    INTEGER, INTENT(IN) :: J(2)
    REAL(KIND=DWP), INTENT(OUT) :: U(2,2), Z(2,2)
  END SUBROUTINE DTRANSFA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE DTRANSF
