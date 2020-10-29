MODULE KTRANSF
  USE PARAMS
  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE KHSVD2T(A, U, Z)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2)
    REAL(KIND=DWP), INTENT(OUT) :: U(2,2), Z(2,2)

    REAL(KIND=DWP) :: X, Y, T2, TU, CU, TZ, CZ

    X = MAX(A(1,2) / A(1,1), D_ZERO)
    Y = MAX(A(2,2) / A(1,1), D_ZERO)
    T2 = MIN(MAX((SCALE(MIN(X, Y), 1) * MAX(X, Y)) / ((X - Y) * (X + Y) + D_ONE), D_ZERO), TWOF)
    TU = T2 / (D_ONE + SQRT(T2 * T2 + D_ONE))
    CU = SQRT(TU * TU + D_ONE) ! D_ONE /
    TZ = Y * TU + X
    CZ = SQRT(TZ * TZ + D_ONE) ! D_ONE /

    A(1,1) = (CZ / CU) * A(1,1)
    A(2,2) = (CU / CZ) * A(2,2)

    U(1,1) =  CU
    U(2,1) = -TU
    U(1,2) =  TU
    U(2,2) =  CU

    Z(1,1) =  CZ
    Z(2,1) =  TZ
    Z(1,2) = -TZ
    Z(2,2) =  CZ
  END SUBROUTINE KHSVD2T

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE KHSVD2U(A, U, Z, INFO)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2)
    REAL(KIND=DWP), INTENT(OUT) :: U(2,2), Z(2,2)
    INTEGER, INTENT(OUT) :: INFO

    REAL(KIND=DWP) :: X, Y, T2, TU, CU, TZ, CZ

    X = MAX(A(1,2) / A(1,1), D_ZERO)
    Y = MAX(A(2,2) / A(1,1), D_ZERO)
    T2 = -MIN(MAX((SCALE(MIN(X, Y), 1) * MAX(X, Y)) / ((Y - X) * (Y + X) + D_ONE), T2, D_ZERO), TWOF)
    TU = T2 / (D_ONE + SQRT(T2 * T2 + D_ONE))
    TZ = -(Y * TU + X)
    IF (ABS(TZ) .LT. D_ONE) THEN
       CU = SQRT(TU * TU + D_ONE) ! D_ONE /
       CZ = SQRT(D_ONE - TZ * TZ) ! D_ONE /

       A(1,1) = (CZ / CU) * A(1,1)
       A(2,2) = (CU / CZ) * A(2,2)

       U(1,1) =  CU
       U(2,1) = -TU
       U(1,2) =  TU
       U(2,2) =  CU

       Z(1,1) =  CZ
       Z(2,1) =  TZ
       Z(1,2) =  TZ
       Z(2,2) =  CZ
    ELSE ! |tanh| >= 1
       INFO = -10
    END IF
  END SUBROUTINE KHSVD2U

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE KHSVD2L(A, U, Z, INFO)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2)
    REAL(KIND=DWP), INTENT(OUT) :: U(2,2), Z(2,2)
    INTEGER, INTENT(OUT) :: INFO

    REAL(KIND=DWP) :: X, Y, T2, TU, CU, TZ, CZ

    X = MAX(A(2,1) / A(2,2), D_ZERO)
    Y = MAX(A(1,1) / A(2,2), D_ZERO)
    T2 = MIN(MAX((SCALE(MIN(X, Y), 1) * MAX(X, Y)) / ((Y - X) * (Y + X) + D_ONE), D_ZERO), TWOF)
    TU = T2 / (D_ONE + SQRT(T2 * T2 + D_ONE))
    TZ = Y * TU - X
    IF (ABS(TZ) .LT. D_ONE) THEN
       CU = SQRT(TU * TU + D_ONE) ! D_ONE /
       CZ = SQRT(D_ONE - TZ * TZ) ! D_ONE /

       A(1,1) = (CU / CZ) * A(1,1)
       A(2,2) = (CZ / CU) * A(2,2)

       U(1,1) =  CU
       U(2,1) = -TU
       U(1,2) =  TU
       U(2,2) =  CU

       Z(1,1) =  CZ
       Z(2,1) =  TZ
       Z(1,2) =  TZ
       Z(2,2) =  CZ
    ELSE ! |tanh| >= 1
       INFO = -11
    END IF
  END SUBROUTINE KHSVD2L

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE KHSVD2(A, J, U, Z, INFO)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2)
    REAL(KIND=DWP), INTENT(OUT) :: U(2,2), Z(2,2)
    INTEGER, INTENT(IN) :: J(2)
    INTEGER, INTENT(INOUT) :: INFO

#ifdef NDEBUG
    IF (A(2,1) .EQ. D_ZERO) THEN ! upper triangular
#else
    IF ((INFO .NE. 0) .AND. (INFO .NE. 1)) THEN ! error
       INFO = -1
    ELSE IF (ABS(J(1)) .NE. 1) THEN ! error
       INFO = -2
    ELSE IF (ABS(J(2)) .NE. 1) THEN ! error
       INFO = -3
    ELSE IF (.NOT. (A(1,1) .GE. D_ZERO)) THEN ! error
       INFO = -4
    ELSE IF (.NOT. (A(2,1) .GE. D_ZERO)) THEN ! error
       INFO = -5
    ELSE IF (.NOT. (A(1,2) .GE. D_ZERO)) THEN ! error
       INFO = -6
    ELSE IF (.NOT. (A(2,2) .GE. D_ZERO)) THEN ! error
       INFO = -7
    ELSE IF (A(2,1) .EQ. D_ZERO) THEN ! upper triangular
#endif
       IF (J(1) .EQ. J(2)) THEN ! trigonometric
          IF (INFO .EQ. 1) CALL KHSVD2T(A, U, Z)
       ELSE ! hyperbolic
          IF (INFO .EQ. 1) CALL KHSVD2U(A, U, Z, INFO)
       END IF
    ELSE IF (A(1,2) .EQ. D_ZERO) THEN ! hyperbolic lower triangular
#ifndef NDEBUG
       IF (J(1) .EQ. J(2)) THEN ! error
          INFO = -8
       ELSE ! non-diagonal
#endif
          IF (INFO .EQ. 1) CALL KHSVD2L(A, U, Z, INFO)
#ifndef NDEBUG
       END IF
    ELSE ! A not triangular => error
       INFO = -9
#endif
    END IF
  END SUBROUTINE KHSVD2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE KTRANSF
