MODULE DTRANSF
  USE PARAMS
  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE DTRANSFC(A11, A21, A12, A22, J1, J2, CU, SU, CZ, SZ)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: A11, A21, A12, A22
    INTEGER, INTENT(IN) :: J1, J2
    REAL(KIND=DWP), INTENT(OUT) :: CU, SU, CZ, SZ

    CU = D_ONE
    SU = D_ZERO
    CZ = D_ONE
    SZ = D_ZERO
  END SUBROUTINE DTRANSFC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE DTRANSFA(A11, A21, A12, A22, J1, J2, CU, SU, CZ, SZ, LEFT, A, U, Z, LI, RI, NC)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: A11, A21, A12, A22
    INTEGER, INTENT(IN) :: J1, J2
    REAL(KIND=DWP), INTENT(OUT) :: CU, CZ, SU, SZ, A(2,2), U(2,2), Z(2,2)
    LOGICAL, INTENT(IN) :: LEFT
    LOGICAL, INTENT(OUT) :: LI, RI, NC ! no change

    REAL(KIND=DWP) :: UA(2,2), AZ(2,2)
    EQUIVALENCE (UA, AZ)

    CALL DTRANSFC(A11, A21, A12, A22, J1, J2, CU, SU, CZ, SZ)

    A(1,1) = A11
    A(2,1) = A21
    A(1,2) = A12
    A(2,2) = A22

    U(1,1) =  CU
    U(2,1) = -SU
    U(1,2) =  SU
    U(2,2) =  CU
    LI = ((CU .EQ. D_ONE) .AND. (SU .EQ. D_ZERO))

    Z(1,1) =  CZ
    Z(2,1) = -SZ
    Z(1,2) =  SZ
    Z(2,2) =  CZ
    RI = ((CZ .EQ. D_ONE) .AND. (SZ .EQ. D_ZERO))

    IF (LEFT) THEN
       IF (LI) THEN
          UA = A
       ELSE ! U .NE. I
          UA = MATMUL(U, A)
       END IF
       IF (RI) THEN
          A = UA
       ELSE ! R .NE. I
          A = MATMUL(UA, Z)
       END IF
    ELSE ! RIGHT
       IF (RI) THEN
          AZ = A
       ELSE ! R .NE. I
          AZ = MATMUL(A, Z)
       END IF
       IF (LI) THEN
          A = AZ
       ELSE ! U .NE. I
          A = MATMUL(U, AZ)
       END IF
    END IF

    NC = ((A11 .EQ. A(1,1)) .AND. (A21 .EQ. A(2,1)) .AND. (A12 .EQ. A(1,2)) .AND. (A22 .EQ. A(2,2)))
  END SUBROUTINE DTRANSFA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE DTRANSF
