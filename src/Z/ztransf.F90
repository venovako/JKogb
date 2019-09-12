MODULE ZTRANSF
  USE PARAMS
  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE ZTRANSFA(A, J, U, Z)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(2,2)
    INTEGER, INTENT(IN) :: J(2)
    COMPLEX(KIND=DWP), INTENT(OUT) :: U(2,2), Z(2,2)
  END SUBROUTINE ZTRANSFA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE ZTRANSF
