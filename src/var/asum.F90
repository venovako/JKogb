!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE REAL(KIND=DWP) FUNCTION DASUM4(X1, X2, X3, X4)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: X1, X2, X3, X4

    DASUM4 = X1 * X1
    DASUM4 = FMAD(X2, X2, DASUM4)
    DASUM4 = FMAD(X3, X3, DASUM4)
    DASUM4 = FMAD(X4, X4, DASUM4)
  END FUNCTION DASUM4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE REAL(KIND=DWP) FUNCTION DASUM2(X1, X2)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: X1, X2

    REAL(KIND=DWP) :: X(2)

    X(1) = ABS(X1)
    X(2) = ABS(X2)

    IF (X(1) .GE. X(2)) THEN
       IF (X(1) .EQ. D_ZERO) THEN
          DASUM2 = D_ZERO
       ELSE ! X(1) .GT. D_ZERO
          DASUM2 = FMAD(X(2), (X(2) / X(1)), X(1)) * X(1)
       END IF
    ELSE ! X(1) .LT. X(2)
       DASUM2 = FMAD(X(1), (X(1) / X(2)), X(2)) * X(2)
    END IF
  END FUNCTION DASUM2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
