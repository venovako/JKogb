MODULE DTRANSF
  USE PARAMS
  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE UT(U)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(INOUT) :: U(2,2)

    REAL(KIND=DWP) :: U21

    U21 = U(2,1)
    U(2,1) = U(1,2)
    U(1,2) = U21
  END SUBROUTINE UT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE CA(B, N, X, Y, LDA)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: B(2,2)
    INTEGER, INTENT(IN) :: N, LDA
    REAL(KIND=DWP), INTENT(INOUT) :: X(*), Y(*)

    REAL(KIND=DWP) :: XX, YY
    INTEGER :: I, J

    I = 1
    !DIR$ VECTOR ALWAYS
    DO J = 1, N
       !DIR$ FMA
       XX = X(I) + B(1,2) * Y(I)
       !DIR$ FMA
       YY = B(2,1) * X(I) + Y(I)
       X(I) = B(1,1) * XX
       Y(I) = B(2,2) * YY
       I = I + LDA
    END DO
  END SUBROUTINE CA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE BA(B, N, X, Y, LDA)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: B(2,2)
    INTEGER, INTENT(IN) :: N, LDA
    REAL(KIND=DWP), INTENT(INOUT) :: X(*), Y(*)

    REAL(KIND=DWP) :: R1, R2, XX, YY
    INTEGER :: I, J

    I = 1

    ! DO J = 1, N
    !    XX = B(1,1) * X(I) + B(1,2) * Y(I)
    !    YY = B(2,1) * X(I) + B(2,2) * Y(I)
    !    X(I) = XX
    !    Y(I) = YY
    !    I = I + LDA
    ! END DO

    IF (ABS(B(1,1)) .GE. ABS(B(1,2))) THEN
       R1 = B(1,2) / B(1,1)
       IF (ABS(B(2,2)) .GE. ABS(B(2,1))) THEN
          R2 = B(2,1) / B(2,2)
          !DIR$ VECTOR ALWAYS
          DO J = 1, N
             !DIR$ FMA
             XX = X(I) + R1 * Y(I)
             !DIR$ FMA
             YY = R2 * X(I) + Y(I)
             X(I) = XX * B(1,1)
             Y(I) = YY * B(2,2)
             I = I + LDA
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(2,1))
          R2 = B(2,2) / B(2,1)
          !DIR$ VECTOR ALWAYS
          DO J = 1, N
             !DIR$ FMA
             XX = X(I) + R1 * Y(I)
             !DIR$ FMA
             YY = X(I) + R2 * Y(I)
             X(I) = XX * B(1,1)
             Y(I) = YY * B(2,1)
             I = I + LDA
          END DO
       END IF
    ELSE ! ABS(B(1,1)) .LT. ABS(B(1,2))
       R1 = B(1,1) / B(1,2)
       IF (ABS(B(2,2)) .GE. ABS(B(2,1))) THEN
          R2 = B(2,1) / B(2,2)
          !DIR$ VECTOR ALWAYS
          DO J = 1, N
             !DIR$ FMA
             XX = R1 * X(I) + Y(I)
             !DIR$ FMA
             YY = R2 * X(I) + Y(I)
             X(I) = XX * B(1,2)
             Y(I) = YY * B(2,2)
             I = I + LDA
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(2,1))
          R2 = B(2,2) / B(2,1)
          !DIR$ VECTOR ALWAYS
          DO J = 1, N
             !DIR$ FMA
             XX = R1 * X(I) + Y(I)
             !DIR$ FMA
             YY = X(I) + R2 * Y(I)
             X(I) = XX * B(1,2)
             Y(I) = YY * B(2,1)
             I = I + LDA
          END DO
       END IF
    END IF
  END SUBROUTINE BA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE AC(B, M, X, Y)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: B(2,2)
    INTEGER, INTENT(IN) :: M
    REAL(KIND=DWP), INTENT(INOUT) :: X(M), Y(M)

    REAL(KIND=DWP) :: XX, YY
    INTEGER :: I

    !DIR$ VECTOR ALWAYS
    DO I = 1, M
       !DIR$ FMA
       XX = X(I) + Y(I) * B(2,1)
       !DIR$ FMA
       YY = X(I) * B(1,2) + Y(I)
       X(I) = XX * B(1,1)
       Y(I) = YY * B(2,2)
    END DO
  END SUBROUTINE AC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE AB(B, M, X, Y)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: B(2,2)
    INTEGER, INTENT(IN) :: M
    REAL(KIND=DWP), INTENT(INOUT) :: X(M), Y(M)

    REAL(KIND=DWP) :: R1, R2, XX, YY
    INTEGER :: I

    ! DO I = 1, M
    !    XX = X(I) * B(1,1) + Y(I) * B(2,1)
    !    YY = X(I) * B(1,2) + Y(I) * B(2,2)
    !    X(I) = XX
    !    Y(I) = YY
    ! END DO

    IF (ABS(B(1,1)) .GE. ABS(B(2,1))) THEN
       R1 = B(2,1) / B(1,1)
       IF (ABS(B(2,2)) .GE. ABS(B(1,2))) THEN
          R2 = B(1,2) / B(2,2)
          !DIR$ VECTOR ALWAYS
          DO I = 1, M
             !DIR$ FMA
             XX = X(I) + Y(I) * R1
             !DIR$ FMA
             YY = X(I) * R2 + Y(I)
             X(I) = XX * B(1,1)
             Y(I) = YY * B(2,2)
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(1,2))
          R2 = B(2,2) / B(1,2)
          !DIR$ VECTOR ALWAYS
          DO I = 1, M
             !DIR$ FMA
             XX = X(I) + Y(I) * R1
             !DIR$ FMA
             YY = X(I) + Y(I) * R2
             X(I) = XX * B(1,1)
             Y(I) = YY * B(1,2)
          END DO
       END IF
    ELSE ! ABS(B(1,1)) .LT. ABS(B(2,1))
       R1 = B(1,1) / B(2,1)
       IF (ABS(B(2,2)) .GE. ABS(B(1,2))) THEN
          R2 = B(1,2) / B(2,2)
          !DIR$ VECTOR ALWAYS
          DO I = 1, M
             !DIR$ FMA
             XX = X(I) * R1 + Y(I)
             !DIR$ FMA
             YY = X(I) * R2 + Y(I)
             X(I) = XX * B(2,1)
             Y(I) = YY * B(2,2)
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(1,2))
          R2 = B(2,2) / B(1,2)
          !DIR$ VECTOR ALWAYS
          DO I = 1, M
             !DIR$ FMA
             XX = X(I) * R1 + Y(I)
             !DIR$ FMA
             YY = X(I) + Y(I) * R2
             X(I) = XX * B(2,1)
             Y(I) = YY * B(1,2)
          END DO
       END IF
    END IF
  END SUBROUTINE AB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE REAL(KIND=DWP) FUNCTION ABODND(B, M, X, Y, P, Q)
    ! skipping the diagonal elements X(P) and Y(Q),
    ! D = SUM(oldX(I)**2-newX(I)**2 + oldY(I)**2-newY(I)**2)
    !   = ||oldX oldY||_F^2 - ||newX newY||_F^2
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: B(2,2), X(M), Y(M)
    INTEGER, INTENT(IN) :: M, P, Q

    REAL(KIND=DWP) :: R1, R2, XX, YY
    INTEGER :: I

    IF ((M .LT. 2) .OR. (P .LE. 0) .OR. (P .GT. M) .OR. (Q .LE. P) .OR. (Q .GT. M)) THEN
       ABODND = D_MZERO
       RETURN
    END IF
    ABODND = Y(P) * Y(P) + X(Q) * X(Q)

    IF (ABS(B(1,1)) .GE. ABS(B(2,1))) THEN
       R1 = B(2,1) / B(1,1)
       IF (ABS(B(2,2)) .GE. ABS(B(1,2))) THEN
          R2 = B(1,2) / B(2,2)
          !DIR$ VECTOR ALWAYS
          DO I = 1, P-1
             !DIR$ FMA
             XX = (X(I) + Y(I) * R1) * B(1,1)
             !DIR$ FMA
             YY = (X(I) * R2 + Y(I)) * B(2,2)
             !DIR$ FMA
             ABODND = ABODND + (X(I) - XX) * (X(I) + XX)
             !DIR$ FMA
             ABODND = ABODND + (Y(I) - YY) * (Y(I) + YY)
          END DO
          !DIR$ VECTOR ALWAYS
          DO I = P+1, Q-1
             !DIR$ FMA
             XX = (X(I) + Y(I) * R1) * B(1,1)
             !DIR$ FMA
             YY = (X(I) * R2 + Y(I)) * B(2,2)
             !DIR$ FMA
             ABODND = ABODND + (X(I) - XX) * (X(I) + XX)
             !DIR$ FMA
             ABODND = ABODND + (Y(I) - YY) * (Y(I) + YY)
          END DO
          !DIR$ VECTOR ALWAYS
          DO I = Q+1, M
             !DIR$ FMA
             XX = (X(I) + Y(I) * R1) * B(1,1)
             !DIR$ FMA
             YY = (X(I) * R2 + Y(I)) * B(2,2)
             !DIR$ FMA
             ABODND = ABODND + (X(I) - XX) * (X(I) + XX)
             !DIR$ FMA
             ABODND = ABODND + (Y(I) - YY) * (Y(I) + YY)
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(1,2))
          R2 = B(2,2) / B(1,2)
          !DIR$ VECTOR ALWAYS
          DO I = 1, P-1
             !DIR$ FMA
             XX = (X(I) + Y(I) * R1) * B(1,1)
             !DIR$ FMA
             YY = (X(I) + Y(I) * R2) * B(1,2)
             !DIR$ FMA
             ABODND = ABODND + (X(I) - XX) * (X(I) + XX)
             !DIR$ FMA
             ABODND = ABODND + (Y(I) - YY) * (Y(I) + YY)
          END DO
          !DIR$ VECTOR ALWAYS
          DO I = P+1, Q-1
             !DIR$ FMA
             XX = (X(I) + Y(I) * R1) * B(1,1)
             !DIR$ FMA
             YY = (X(I) + Y(I) * R2) * B(1,2)
             !DIR$ FMA
             ABODND = ABODND + (X(I) - XX) * (X(I) + XX)
             !DIR$ FMA
             ABODND = ABODND + (Y(I) - YY) * (Y(I) + YY)
          END DO
          !DIR$ VECTOR ALWAYS
          DO I = Q+1, M
             !DIR$ FMA
             XX = (X(I) + Y(I) * R1) * B(1,1)
             !DIR$ FMA
             YY = (X(I) + Y(I) * R2) * B(1,2)
             !DIR$ FMA
             ABODND = ABODND + (X(I) - XX) * (X(I) + XX)
             !DIR$ FMA
             ABODND = ABODND + (Y(I) - YY) * (Y(I) + YY)
          END DO
       END IF
    ELSE ! ABS(B(1,1)) .LT. ABS(B(2,1))
       R1 = B(1,1) / B(2,1)
       IF (ABS(B(2,2)) .GE. ABS(B(1,2))) THEN
          R2 = B(1,2) / B(2,2)
          !DIR$ VECTOR ALWAYS
          DO I = 1, P-1
             !DIR$ FMA
             XX = (X(I) * R1 + Y(I)) * B(2,1)
             !DIR$ FMA
             YY = (X(I) * R2 + Y(I)) * B(2,2)
             !DIR$ FMA
             ABODND = ABODND + (X(I) - XX) * (X(I) + XX)
             !DIR$ FMA
             ABODND = ABODND + (Y(I) - YY) * (Y(I) + YY)
          END DO
          !DIR$ VECTOR ALWAYS
          DO I = P+1, Q-1
             !DIR$ FMA
             XX = (X(I) * R1 + Y(I)) * B(2,1)
             !DIR$ FMA
             YY = (X(I) * R2 + Y(I)) * B(2,2)
             !DIR$ FMA
             ABODND = ABODND + (X(I) - XX) * (X(I) + XX)
             !DIR$ FMA
             ABODND = ABODND + (Y(I) - YY) * (Y(I) + YY)
          END DO
          !DIR$ VECTOR ALWAYS
          DO I = Q+1, M
             !DIR$ FMA
             XX = (X(I) * R1 + Y(I)) * B(2,1)
             !DIR$ FMA
             YY = (X(I) * R2 + Y(I)) * B(2,2)
             !DIR$ FMA
             ABODND = ABODND + (X(I) - XX) * (X(I) + XX)
             !DIR$ FMA
             ABODND = ABODND + (Y(I) - YY) * (Y(I) + YY)
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(1,2))
          R2 = B(2,2) / B(1,2)
          !DIR$ VECTOR ALWAYS
          DO I = 1, P-1
             !DIR$ FMA
             XX = (X(I) * R1 + Y(I)) * B(2,1)
             !DIR$ FMA
             YY = (X(I) + Y(I) * R2) * B(1,2)
             !DIR$ FMA
             ABODND = ABODND + (X(I) - XX) * (X(I) + XX)
             !DIR$ FMA
             ABODND = ABODND + (Y(I) - YY) * (Y(I) + YY)
          END DO
          !DIR$ VECTOR ALWAYS
          DO I = P+1, Q-1
             !DIR$ FMA
             XX = (X(I) * R1 + Y(I)) * B(2,1)
             !DIR$ FMA
             YY = (X(I) + Y(I) * R2) * B(1,2)
             !DIR$ FMA
             ABODND = ABODND + (X(I) - XX) * (X(I) + XX)
             !DIR$ FMA
             ABODND = ABODND + (Y(I) - YY) * (Y(I) + YY)
          END DO
          !DIR$ VECTOR ALWAYS
          DO I = Q+1, M
             !DIR$ FMA
             XX = (X(I) * R1 + Y(I)) * B(2,1)
             !DIR$ FMA
             YY = (X(I) + Y(I) * R2) * B(1,2)
             !DIR$ FMA
             ABODND = ABODND + (X(I) - XX) * (X(I) + XX)
             !DIR$ FMA
             ABODND = ABODND + (Y(I) - YY) * (Y(I) + YY)
          END DO
       END IF
    END IF
  END FUNCTION ABODND

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE DHSVD2D(H, A, U, Z, INFO)
    ! A diagonal
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: H
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    INFO = 0

    IF (SIGN(D_ONE, A(1,1)) .EQ. D_MONE) THEN
       ! A(1,1) negative
       U(1,1) = -U(1,1)
       U(1,2) = -U(1,2)
       A(1,1) = -A(1,1)
    END IF

    IF (SIGN(D_ONE, A(2,2)) .EQ. D_MONE) THEN
       ! A(2,2) negative
       U(2,1) = -U(2,1)
       U(2,2) = -U(2,2)
       A(2,2) = -A(2,2)
    END IF
  END SUBROUTINE DHSVD2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DHSVD2U(H, A, U, Z, INFO)
    ! A upper triangular, not diagonal
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: H
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    REAL(KIND=DWP) :: V(2,2), W(2,2)
    REAL(KIND=DWP) :: CU, SU, TU ! always trigonometric
    REAL(KIND=DWP) :: CZ, SZ, TZ ! trigonometric or hyperbolic

    REAL(KIND=DWP) :: T2U, X, Y

    INFO = 0

    IF (H) THEN
       IF (A(2,2) .NE. D_ZERO) THEN
          X = A(1,2) / A(1,1) ! .NE. 0
          Y = A(2,2) / A(1,1) ! .GT. 0
          IF (ABS(X) .LE. Y) THEN
             T2U = SCALE(X, 1) * Y
          ELSE ! ABS(X) .GT. Y
             T2U = SCALE(Y, 1) * X
          END IF
          !DIR$ FMA
          T2U = T2U / (D_ONE + (Y - X) * (Y + X))
          IF (.NOT. (ABS(T2U) .LE. HUGE(D_ZERO))) THEN
             TU = SIGN(D_ONE, T2U)
          ELSE ! should always happen
             !DIR$ FMA
             TU = T2U / (D_ONE + SQRT(D_ONE + T2U * T2U))
          END IF
          !DIR$ FMA
          CU = D_ONE / SQRT(D_ONE + TU * TU)
          ! SU = TU * CU
          !DIR$ FMA
          TZ = Y * TU - X
          IF (ABS(TZ) .GE. D_ONE) THEN
             INFO = -6
             RETURN
          END IF
          !DIR$ FMA
          CZ = D_ONE / SQRT(D_ONE - TZ * TZ)
          IF (.NOT. (CZ .LE. HUGE(D_ZERO))) THEN
             INFO = -7
             RETURN
          END IF
          SZ = TZ * CZ
          SZ = SZ * SZ ! SZ^2
          IF ((SZ + D_ONE) .EQ. SZ) THEN
             INFO = -8
             RETURN
          END IF
       ELSE IF (ABS(A(1,2)) .LT. ABS(A(1,1))) THEN
          INFO = 1
          ! TU = D_ZERO
          ! CU = D_ONE
          ! SU = D_ZERO
          TZ = -A(1,2) / A(1,1)
          IF (ABS(TZ) .GE. D_ONE) THEN
             INFO = -9
             RETURN
          END IF
          !DIR$ FMA
          CZ = D_ONE / SQRT(D_ONE - TZ * TZ)
          IF (.NOT. (CZ .LE. HUGE(D_ZERO))) THEN
             INFO = -10
             RETURN
          END IF
          SZ = TZ * CZ
          SZ = SZ * SZ ! SZ^2
          IF ((SZ + D_ONE) .EQ. SZ) THEN
             INFO = -11
             RETURN
          END IF
       ELSE ! (A(2,2) .EQ. 0) .AND. (ABS(A(1,1)) .EQ. ABS(A(1,2)))
          ! |TZ| .GE. 1
          INFO = -12
          RETURN
       END IF
       W(1,1) =  CZ
       W(2,1) =  TZ ! SZ
       W(1,2) =  TZ ! SZ
       W(2,2) =  CZ
    ELSE ! trigonometric
       IF (A(2,2) .NE. D_ZERO) THEN
          X = A(1,2) / A(1,1) ! .NE. 0
          Y = A(2,2) / A(1,1) ! .GT. 0
          IF (ABS(X) .LE. Y) THEN
             T2U = -SCALE(X, 1) * Y
          ELSE ! ABS(X) .GT. Y
             T2U = -SCALE(Y, 1) * X
          END IF
          !DIR$ FMA
          T2U = T2U / (D_ONE + (X - Y) * (X + Y))
          IF (.NOT. (ABS(T2U) .LE. HUGE(D_ZERO))) THEN
             TU = SIGN(D_ONE, T2U)
          ELSE ! should always happen
             !DIR$ FMA
             TU = T2U / (D_ONE + SQRT(D_ONE + T2U * T2U))
          END IF
          !DIR$ FMA
          CU = D_ONE / SQRT(D_ONE + TU * TU)
          ! SU = TU * CU
          !DIR$ FMA
          TZ = Y * TU - X
          !DIR$ FMA
          CZ = D_ONE / SQRT(D_ONE + TZ * TZ)
          ! SZ = TZ * CZ
       ELSE ! A(2,2) .EQ. 0
          INFO = 1
          ! TU = D_ZERO
          ! CU = D_ONE
          ! SU = D_ZERO
          TZ = -A(1,2) / A(1,1)
          !DIR$ FMA
          CZ = D_ONE / SQRT(D_ONE + TZ * TZ)
          ! SZ = TZ * CZ
       END IF
       W(1,1) =  CZ
       W(2,1) = -TZ !-SZ
       W(1,2) =  TZ ! SZ
       W(2,2) =  CZ
    END IF

    IF (INFO .EQ. 0) THEN
       V(1,1) =  CU
       V(2,1) =  TU ! SU
       V(1,2) = -TU !-SU
       V(2,2) =  CU

       CALL CA(V, 2, U(1,1), U(2,1), 2)
       CALL CA(V, 2, A(1,1), A(2,1), 2)
    ELSE ! INFO .NE. 0
       INFO = 0
    END IF

    CALL AC(W, 2, A(1,1), A(1,2))
    CALL AC(W, 2, Z(1,1), Z(1,2))

    CALL DHSVD2D(H, A, U, Z, INFO)
  END SUBROUTINE DHSVD2U

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DHSVD2L(H, A, U, Z, INFO)
    ! A lower triangular, not diagonal
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: H
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    REAL(KIND=DWP) :: V(2,2), W(2,2)
    REAL(KIND=DWP) :: CU, SU, TU ! always trigonometric
    REAL(KIND=DWP) :: CZ, SZ, TZ ! always hyperbolic

    REAL(KIND=DWP) :: T2U, X, Y

    IF (.NOT. H) THEN
       INFO = -HUGE(0)
       RETURN
    END IF
    INFO = 0

    IF (A(1,1) .NE. D_ZERO) THEN
       X = A(1,1) / A(2,2) ! .GT. 0
       Y = A(2,1) / A(2,2) ! .NE. 0
       IF (ABS(Y) .LE. X) THEN
          T2U = -SCALE(Y, 1) * X
       ELSE ! ABS(Y) .GT. X
          T2U = -SCALE(X, 1) * Y
       END IF
       !DIR$ FMA
       T2U = T2U / (D_ONE + (X - Y) * (X + Y))
       IF (.NOT. (ABS(T2U) .LE. HUGE(D_ZERO))) THEN
          TU = SIGN(D_ONE, T2U)
       ELSE ! should always happen
          !DIR$ FMA
          TU = T2U / (D_ONE + SQRT(D_ONE + T2U * T2U))
       END IF
       !DIR$ FMA
       CU = D_ONE / SQRT(D_ONE + TU * TU)
       ! SU = TU * CU
       !DIR$ FMA
       TZ = -(X * TU + Y)
       IF (ABS(TZ) .GE. D_ONE) THEN
          INFO = -13
          RETURN
       END IF
       !DIR$ FMA
       CZ = D_ONE / SQRT(D_ONE - TZ * TZ)
       IF (.NOT. (CZ .LE. HUGE(D_ZERO))) THEN
          INFO = -14
          RETURN
       END IF
       SZ = TZ * CZ
       SZ = SZ * SZ ! SZ^2
       IF ((SZ + D_ONE) .EQ. SZ) THEN
          INFO = -15
          RETURN
       END IF
    ELSE IF (ABS(A(2,1)) .LT. ABS(A(2,2))) THEN
       INFO = 1
       ! TU = D_ZERO
       ! CU = D_ONE
       ! SU = D_ZERO
       TZ = -A(2,1) / A(2,2)
       IF (ABS(TZ) .GE. D_ONE) THEN
          INFO = -16
          RETURN
       END IF
       !DIR$ FMA
       CZ = D_ONE / SQRT(D_ONE - TZ * TZ)
       IF (.NOT. (CZ .LE. HUGE(D_ZERO))) THEN
          INFO = -17
          RETURN
       END IF
       SZ = TZ * CZ
       SZ = SZ * SZ ! SZ^2
       IF ((SZ + D_ONE) .EQ. SZ) THEN
          INFO = -18
          RETURN
       END IF
    ELSE ! (A(1,1) .EQ. 0) .AND. (ABS(A(2,1)) .EQ. ABS(A(2,2)))
       ! |TZ| .GE. 1
       INFO = -19
       RETURN
    END IF

    IF (INFO .EQ. 0) THEN
       V(1,1) =  CU
       V(2,1) =  TU ! SU
       V(1,2) = -TU !-SU
       V(2,2) =  CU

       CALL CA(V, 2, U(1,1), U(2,1), 2)
       CALL CA(V, 2, A(1,1), A(2,1), 2)
    ELSE ! INFO .NE. 0
       INFO = 0
    END IF

    W(1,1) =  CZ
    W(2,1) =  TZ ! SZ
    W(1,2) =  TZ ! SZ
    W(2,2) =  CZ

    CALL AC(W, 2, A(1,1), A(1,2))
    CALL AC(W, 2, Z(1,1), Z(1,2))

    CALL DHSVD2D(H, A, U, Z, INFO)
  END SUBROUTINE DHSVD2L

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DHSVD2G(H, A, U, Z, INFO)
    ! A general
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: H
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    REAL(KIND=DWP) :: C, S, R

    REAL(KIND=DWP), EXTERNAL :: DNRM2
    EXTERNAL :: DLARTG, DROT, DSWAP

    ! column pivoting
    IF (DNRM2(2, A(1,1), 1) .LT. DNRM2(2, A(1,2), 1)) THEN
       ! column swap of A
       CALL DSWAP(2, A(1,1), 1, A(1,2), 1)
       ! column swap of Z
       IF (.NOT. H) CALL DSWAP(2, Z(1,1), 1, Z(1,2), 1)
       ! record the swap
       INFO = 1
    ELSE ! no pivoting
       INFO = 0
    END IF

    ! row sorting
    IF (ABS(A(1,1)) .LT. ABS(A(2,1))) THEN
       ! row swap of U
       CALL DSWAP(2, U(1,1), 2, U(2,1), 2)
       ! row swap of A
       CALL DSWAP(2, A(1,1), 2, A(2,1), 2)
    END IF

    ! QR factorization of A
    CALL DLARTG(A(1,1), A(2,1), C, S, R)
    IF (.NOT. (ABS(R) .LE. HUGE(D_ZERO))) THEN
       INFO = -5
       RETURN
    END IF
    CALL DROT(1, A(1,2), 2, A(2,2), 2, C, S)

    ! premultiply U by Q^T
    CALL DROT(2, U(1,1), 2, U(2,1), 2, C, S)

    ! make diag(A) non-negative
    IF (SIGN(D_ONE, R) .EQ. D_MONE) THEN
       U(1,1) = -U(1,1)
       U(1,2) = -U(1,2)
       A(1,2) = -A(1,2)
       R = -R
    END IF
    A(1,1) = R
    IF (SIGN(D_ONE, A(2,2)) .EQ. D_MONE) THEN
       U(2,1) = -U(2,1)
       U(2,2) = -U(2,2)
       A(2,2) = -A(2,2)
    END IF
    A(2,1) = D_ZERO

    IF (A(1,2) .EQ. D_ZERO) THEN
       ! A diagonal
       IF (H .AND. (INFO .EQ. 1)) THEN
          ! swap the rows of U
          CALL DSWAP(2, U(1,1), 2, U(2,1), 2)
          ! swap the diagonal elements of A
          CALL DSWAP(1, A(1,1), 1, A(2,2), 1)
          ! early exit
          INFO = 0
       END IF
    ELSE IF (H .AND. (INFO .EQ. 1)) THEN
       ! column swap of A
       CALL DSWAP(2, A(1,1), 1, A(1,2), 1)
       ! A upper antitriangular, not antidiagonal (X .NE. 0)
       !     | X R | <- R .NE. 0 the largest
       ! A = | x 0 |    element by magnitude
       ! swap the rows of U
       CALL DSWAP(2, U(1,1), 2, U(2,1), 2)
       ! swap the rows of A
       CALL DSWAP(2, A(1,1), 2, A(2,1), 2)
       ! A lower triangular, not diagonal (X .NE. 0)
       ! A = | x 0 |    R .NE. 0 the largest
       !     | X R | <- element by magnitude
       CALL DHSVD2L(H, A, U, Z, INFO)
    ELSE
       ! A upper triangular, not diagonal (X .NE. 0)
       !     | R X | <- R .NE. 0 the largest
       ! A = | 0 x |    element by magnitude
       CALL DHSVD2U(H, A, U, Z, INFO)
    END IF
  END SUBROUTINE DHSVD2G

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DHSVD2S(H, A, U, Z, INFO)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: H
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    EXTERNAL :: DSWAP

    ! assume A real, non-negative, diagonal
    ! permute A, U, Z in trigonometric case
    IF ((.NOT. H) .AND. (A(1,1) .LT. A(2,2))) THEN
       ! swap the rows of U
       CALL DSWAP(2, U(1,1), 2, U(2,1), 2)
       ! swap the diagonal elements of A
       CALL DSWAP(1, A(1,1), 1, A(2,2), 1)
       ! swap the columns of Z
       CALL DSWAP(2, Z(1,1), 1, Z(1,2), 1)
    END IF

    ! record in INFO if old diag(A) and new diag(A) differ
    IF ((A(1,1) .NE. A(2,1)) .OR. (A(2,2) .NE. A(1,2))) THEN
       INFO = 1
    ELSE ! old diag(A) .EQ. new diag(A)
       INFO = 0
    END IF

    ! explicitly make A diagonal
    A(2,1) = D_ZERO
    A(1,2) = D_ZERO

    ! check if U is identity and record in INFO if it is not
    IF ((U(1,1) .NE. D_ONE) .OR. (U(2,1) .NE. D_ZERO) .OR. (U(1,2) .NE. D_ZERO) .OR. (U(2,2) .NE. D_ONE)) INFO = INFO + 2
    ! check if Z is identity and record in INFO if it is not
    IF ((Z(1,1) .NE. D_ONE) .OR. (Z(2,1) .NE. D_ZERO) .OR. (Z(1,2) .NE. D_ZERO) .OR. (Z(2,2) .NE. D_ONE)) INFO = INFO + 4
  END SUBROUTINE DHSVD2S

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE DSCALEA(A, S)
    IMPLICIT NONE

    REAL(KIND=DWP), PARAMETER :: TOOBIG = SCALE(HUGE(D_ZERO), -1)
    REAL(KIND=DWP), PARAMETER :: TOOSMALL = TINY(D_ZERO)

    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2)
    INTEGER, INTENT(OUT) :: S

    REAL(KIND=DWP) :: AA, AX
    INTEGER :: DS, US, I, J

    DS = 0
    US = 0

    DO J = 1, 2
       DO I = 1, 2
          AA = ABS(A(I,J))
          IF (.NOT. (AA .LE. HUGE(D_ZERO))) THEN
             ! -2 <= S <= -5
             S = -((J - 1) * 2 + I + 1)
             RETURN
          END IF
          ! DS cannot be less than -1, but...
          IF (AA .GE. TOOBIG) DS = MIN(DS, -1)
       END DO
    END DO

    IF (DS .NE. 0) THEN
       DO J = 1, 2
          DO I = 1, 2
             A(I,J) = SCALE(A(I,J), DS)
          END DO
       END DO

       S = DS
    ELSE ! might perform upscaling
       AX = D_ZERO

       DO J = 1, 2
          DO I = 1, 2
             AA = ABS(A(I,J))
             IF (AA .GT. AX) AX = AA
             IF (AA .LT. TOOSMALL) US = MAX(US, (EXPONENT(TOOSMALL) - EXPONENT(AA)))
          END DO
       END DO

       ! how much room there is between AX and TOOBIG
       DS = EXPONENT(TOOBIG) - EXPONENT(AX)
       IF (DS .LE. 0) THEN
          ! DS cannot be less than 0, but...
          US = 0
       ELSE IF (US .GE. DS) THEN
          ! now DS > 0
          IF (SCALE(AX, DS) .GE. TOOBIG) THEN
             ! can only be .LE., but it must not be .EQ.
             US = DS - 1
          ELSE ! cannot reach TOOBIG
             US = DS
          END IF
       ELSE ! US .LT. DS, so US is fine
          CONTINUE
       END IF

       IF (US .NE. 0) THEN
          DO J = 1, 2
             DO I = 1, 2
                A(I,J) = SCALE(A(I,J), US)
             END DO
          END DO
       END IF

       S = US
    END IF
  END SUBROUTINE DSCALEA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DHSVD2(H, A, U, Z, INFO)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: H
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2)
    REAL(KIND=DWP), INTENT(OUT) :: U(2,2), Z(2,2)
    INTEGER, INTENT(OUT) :: INFO

    REAL(KIND=DWP) :: W(2,2)
    INTEGER :: S

    W(2,1) = A(1,1)
    W(1,2) = A(2,2)

    CALL DSCALEA(A, S)
    IF (S .LT. -1) THEN
       ! A has NaNs and/or infinities
       INFO = S + 1
       RETURN
    END IF
    INFO = 0

    ! U = I
    U(1,1) = D_ONE
    U(2,1) = D_ZERO
    U(1,2) = D_ZERO
    U(2,2) = D_ONE

    ! Z = I
    Z(1,1) = D_ONE
    Z(2,1) = D_ZERO
    Z(1,2) = D_ZERO
    Z(2,2) = D_ONE

    IF ((A(2,1) .EQ. D_ZERO) .AND. (A(1,2) .EQ. D_ZERO)) THEN
       ! A diagonal
       CALL DHSVD2D(H, A, U, Z, INFO)
    ELSE
       ! A general
       CALL DHSVD2G(H, A, U, Z, INFO)
    END IF
    IF (INFO .LT. 0) RETURN

    ! scale back if necessary
    IF (S .NE. 0) THEN
       A(1,1) = SCALE(A(1,1), -S)
       A(2,2) = SCALE(A(2,2), -S)
    END IF
    A(2,1) = W(2,1)
    A(1,2) = W(1,2)

    CALL DHSVD2S(H, A, U, Z, INFO)
  END SUBROUTINE DHSVD2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE DTRANSF
