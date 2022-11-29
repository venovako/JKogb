MODULE ktransf
  USE utils
#ifndef FMAD
#define FMAD(a,b,c) ((a)*(b)+(c))
#endif
#ifndef MAXD
#define MAXD MAX
#endif
  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE KHSVD2T(A, G, U, Z, S)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: A(2,2)
    REAL(KIND=DWP), INTENT(OUT) :: G(2), U(2,2), Z(2,2)
    INTEGER, INTENT(IN) :: S

    REAL(KIND=DWP) :: X, Y, T2, TU, CU, TZ, CZ
    INTEGER :: D, E

    X = MAXD(A(1,2) / A(1,1), D_ZERO)
    Y = MAXD(A(2,2) / A(1,1), D_ZERO)
    T2 = MIN(MAXD((SCALE(MIN(X, Y), 1) * MAX(X, Y)) / FMAD((X - Y), (X + Y), D_ONE), D_ZERO), TWOF)
    TU = T2 / (D_ONE + SQRT(FMAD(T2, T2, D_ONE)))
    CU = SQRT(FMAD(TU, TU, D_ONE)) ! D_ONE /
    TZ = FMAD(Y, TU, X)
    CZ = SQRT(FMAD(TZ, TZ, D_ONE)) ! D_ONE /

    Y = CZ / CU
    X = Y
    D = EXPONENT(X) - ET
    IF (S .GT. D) THEN
       X = SCALE(X, -D)
       E = S - D
    ELSE ! full downscaling
       X = SCALE(X, -S)
       E = 0
    END IF

    G(1) = SCALE(A(1,1), -E) * X
    G(2) = SCALE(A(2,2) / Y, -S)

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

  PURE SUBROUTINE KHSVD2U(A, G, U, Z, S, TH1, INFO)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: A(2,2), TH1
    REAL(KIND=DWP), INTENT(OUT) :: G(2), U(2,2), Z(2,2)
    INTEGER, INTENT(IN) :: S
    INTEGER, INTENT(OUT) :: INFO

    REAL(KIND=DWP) :: X, Y, T2, TU, CU, TZ, CZ
    INTEGER :: D, E

    X = MAXD(A(1,2) / A(1,1), D_ZERO)
    Y = MAXD(A(2,2) / A(1,1), D_ZERO)
    T2 = -MIN(MAXD((SCALE(MIN(X, Y), 1) * MAX(X, Y)) / FMAD((Y - X), (Y + X), D_ONE), D_ZERO), TWOF)
    TU = T2 / (D_ONE + SQRT(FMAD(T2, T2, D_ONE)))
    TZ = -FMAD(Y, TU, X)

    T2 = ABS(TZ)
    IF (T2 .GE. D_ONE) THEN
       INFO = -10
       RETURN
    END IF

    IF (T2 .GE. TH1) THEN
       TZ = SIGN(TH1FIX, TZ)
       INFO = INFO + 8
    END IF
    CU = SQRT(FMAD( TU, TU, D_ONE)) ! D_ONE /
    CZ = SQRT(FMAD(-TZ, TZ, D_ONE)) ! D_ONE /

    IF (T2 .LT. TH1) THEN
       Y = CU / CZ
       X = Y
       D = EXPONENT(X) - ET
       IF (S .GT. D) THEN
          X = SCALE(X, -D)
          E = S - D
       ELSE ! full downscaling
          X = SCALE(X, -S)
          E = 0
       END IF
       G(1) = SCALE(A(1,1) / Y, -S)
       G(2) = SCALE(A(2,2), -E) * X
    END IF

    U(1,1) =  CU
    U(2,1) = -TU
    U(1,2) =  TU
    U(2,2) =  CU

    Z(1,1) =  CZ
    Z(2,1) =  TZ
    Z(1,2) =  TZ
    Z(2,2) =  CZ
  END SUBROUTINE KHSVD2U

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE KHSVD2L(A, G, U, Z, S, TH1, INFO)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: A(2,2), TH1
    REAL(KIND=DWP), INTENT(OUT) :: G(2), U(2,2), Z(2,2)
    INTEGER, INTENT(IN) :: S
    INTEGER, INTENT(OUT) :: INFO

    REAL(KIND=DWP) :: X, Y, T2, TU, CU, TZ, CZ
    INTEGER :: D, E

    X = MAXD(A(2,1) / A(2,2), D_ZERO)
    Y = MAXD(A(1,1) / A(2,2), D_ZERO)
    T2 = MIN(MAXD((SCALE(MIN(X, Y), 1) * MAX(X, Y)) / FMAD((Y - X), (Y + X), D_ONE), D_ZERO), TWOF)
    TU = T2 / (D_ONE + SQRT(FMAD(T2, T2, D_ONE)))
    TZ = FMAD(Y, TU, -X)

    T2 = ABS(TZ)
    IF (T2 .GE. D_ONE) THEN
       INFO = -11
       RETURN
    END IF

    IF (T2 .GE. TH1) THEN
       TZ = SIGN(TH1FIX, TZ)
       INFO = INFO + 8
    END IF
    CU = SQRT(FMAD( TU, TU, D_ONE)) ! D_ONE /
    CZ = SQRT(FMAD(-TZ, TZ, D_ONE)) ! D_ONE /

    IF (T2 .LT. TH1) THEN
       Y = CU / CZ
       X = Y
       D = EXPONENT(X) - ET
       IF (S .GT. D) THEN
          X = SCALE(X, -D)
          E = S - D
       ELSE ! full downscaling
          X = SCALE(X, -S)
          E = 0
       END IF
       G(1) = SCALE(A(1,1), -E) * X
       G(2) = SCALE(A(2,2) / Y, -S)
    END IF

    U(1,1) =  CU
    U(2,1) = -TU
    U(1,2) =  TU
    U(2,2) =  CU

    Z(1,1) =  CZ
    Z(2,1) =  TZ
    Z(1,2) =  TZ
    Z(2,2) =  CZ
  END SUBROUTINE KHSVD2L

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE KHSVD2(A, J, G, U, Z, S, TH1, INFO)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: A(2,2), TH1
    REAL(KIND=DWP), INTENT(OUT) :: G(2), U(2,2), Z(2,2)
    INTEGER, INTENT(IN) :: J(2), S
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
          IF (INFO .EQ. 1) CALL KHSVD2T(A, G, U, Z, S)
       ELSE ! hyperbolic
          IF (INFO .EQ. 1) CALL KHSVD2U(A, G, U, Z, S, TH1, INFO)
       END IF
    ELSE IF (A(1,2) .EQ. D_ZERO) THEN ! hyperbolic lower triangular
#ifndef NDEBUG
       IF (J(1) .EQ. J(2)) THEN ! error
          INFO = -8
       ELSE ! non-diagonal
#endif
          IF (INFO .EQ. 1) CALL KHSVD2L(A, G, U, Z, S, TH1, INFO)
#ifndef NDEBUG
       END IF
    ELSE ! A not triangular => error
       INFO = -9
#endif
    END IF
  END SUBROUTINE KHSVD2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE ktransf
