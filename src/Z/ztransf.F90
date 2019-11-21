MODULE ZTRANSF
  USE PARAMS
  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE COMPLEX(KIND=DWP) FUNCTION ZDFMA(A, B, C)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: A
    COMPLEX(KIND=DWP), INTENT(IN) :: B, C

    REAL(KIND=DWP) :: R, I

    R = A * REAL(B) + REAL(C)
    I = A * AIMAG(B) + AIMAG(C)

    ZDFMA = CMPLX(R, I, DWP)
  END FUNCTION ZDFMA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE COMPLEX(KIND=DWP) FUNCTION ZJFMA(A, B, C)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: A ! imaginary
    COMPLEX(KIND=DWP), INTENT(IN) :: B, C

    REAL(KIND=DWP) :: R, I

    R = -A * AIMAG(B) + REAL(C)
    I = A * REAL(B) + AIMAG(C)

    ZJFMA = CMPLX(R, I, DWP)
  END FUNCTION ZJFMA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! after cuCfma() from CUDA's cuComplex.h

  PURE COMPLEX(KIND=DWP) FUNCTION ZZFMA(A, B, C)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(IN) :: A, B, C

    REAL(KIND=DWP) :: F, R, I

    F = REAL(C) - AIMAG(A) * AIMAG(B)
    R = REAL(A) * REAL(B) + F
    F = AIMAG(A) * REAL(B) + AIMAG(C)
    I = REAL(A) * AIMAG(B) + F

    ZZFMA = CMPLX(R, I, DWP)
  END FUNCTION ZZFMA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE COMPLEX(KIND=DWP) FUNCTION ZFMA(A, B, C)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(IN) :: A, B, C

    IF (AIMAG(A) .EQ. D_ZERO) THEN
       ZFMA = ZDFMA(REAL(A), B, C)
    ELSE IF (REAL(A) .EQ. D_ZERO) THEN
       ZFMA = ZJFMA(AIMAG(A), B, C)
    ELSE ! A complex
       ZFMA = ZZFMA(A, B, C)
    END IF
  END FUNCTION ZFMA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE COMPLEX(KIND=DWP) FUNCTION ZZMUL(A, B)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(IN) :: A, B

    REAL(KIND=DWP) :: F, R, I

    F = AIMAG(A) * AIMAG(B)
    R = REAL(A) * REAL(B) - F
    F = AIMAG(A) * REAL(B)
    I = REAL(A) * AIMAG(B) + F

    ZZMUL = CMPLX(R, I, DWP)
  END FUNCTION ZZMUL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! assume that B is not real and |B| .NE. 0
  PURE COMPLEX(KIND=DWP) FUNCTION ZZDIV(A, B)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(IN) :: A, B
#ifdef NDEBUG
    ZZDIV = A / B
#else
    REAL(KIND=DWP) :: F

    F = ABS(B)
    ZZDIV = ZZMUL(CONJG(B / F), A) / F
#endif
  END FUNCTION ZZDIV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE UH2(U)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(INOUT) :: U(2,2)

    COMPLEX(KIND=DWP) :: U21

    U(1,1) = CONJG(U(1,1))
    U21 = U(2,1)
    U(2,1) = CONJG(U(1,2))
    U(1,2) = CONJG(U21)
    U(2,2) = CONJG(U(2,2))
  END SUBROUTINE UH2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE JZHJ2(Z, J)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(INOUT) :: Z(2,2)
    INTEGER, INTENT(IN) :: J(2)

    CALL UH2(Z)
    IF (J(1) .EQ. -1) THEN
       Z(2,1) = -Z(2,1)
       Z(1,2) = -Z(1,2)
    END IF
    IF (J(2) .EQ. -1) THEN
       Z(2,1) = -Z(2,1)
       Z(1,2) = -Z(1,2)
    END IF
  END SUBROUTINE JZHJ2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE CA(B, N, X, Y, LDA)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(IN) :: B(2,2)
    INTEGER, INTENT(IN) :: N, LDA
    COMPLEX(KIND=DWP), INTENT(INOUT) :: X(*), Y(*)

    COMPLEX(KIND=DWP) :: XX, YY
    INTEGER :: I, J

    I = 1

    ! DO J = 1, N
    !    XX = X(I) + B(1,2) * Y(I)
    !    YY = B(2,1) * X(I) + Y(I)
    !    X(I) = B(1,1) * XX
    !    Y(I) = B(2,2) * YY
    !    I = I + LDA
    ! END DO

    DO J = 1, N
       XX = ZFMA(B(1,2), Y(I), X(I))
       YY = ZFMA(B(2,1), X(I), Y(I))
       X(I) = REAL(B(1,1)) * XX
       Y(I) = REAL(B(2,2)) * YY
       I = I + LDA
    END DO
  END SUBROUTINE CA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE BA(B, N, X, Y, LDA)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(IN) :: B(2,2)
    INTEGER, INTENT(IN) :: N, LDA
    COMPLEX(KIND=DWP), INTENT(INOUT) :: X(*), Y(*)

    COMPLEX(KIND=DWP) :: R1, R2, XX, YY
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
       IF (AIMAG(B(1,1)) .EQ. D_ZERO) THEN
          R1 = B(1,2) / REAL(B(1,1))
       ELSE ! B(1,1) complex
          R1 = ZZDIV(B(1,2), B(1,1))
       END IF
       IF (ABS(B(2,2)) .GE. ABS(B(2,1))) THEN
          IF (AIMAG(B(2,2)) .EQ. D_ZERO) THEN
             R2 = B(2,1) / REAL(B(2,2))
          ELSE ! B(2,2) complex
             R2 = ZZDIV(B(2,1), B(2,2))
          END IF
          DO J = 1, N
             XX = ZZFMA(R1, Y(I), X(I)) !X(I) + R1 * Y(I)
             YY = ZZFMA(R2, X(I), Y(I)) !R2 * X(I) + Y(I)
             X(I) = ZZMUL(XX, B(1,1))
             Y(I) = ZZMUL(YY, B(2,2))
             I = I + LDA
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(2,1))
          IF (AIMAG(B(2,1)) .EQ. D_ZERO) THEN
             R2 = B(2,2) / REAL(B(2,1))
          ELSE ! B(2,1) complex
             R2 = ZZDIV(B(2,2), B(2,1))
          END IF
          DO J = 1, N
             XX = ZZFMA(R1, Y(I), X(I)) !X(I) + R1 * Y(I)
             YY = ZZFMA(R1, Y(I), X(I)) !X(I) + R2 * Y(I)
             X(I) = ZZMUL(XX, B(1,1))
             Y(I) = ZZMUL(YY, B(2,1))
             I = I + LDA
          END DO
       END IF
    ELSE ! ABS(B(1,1)) .LT. ABS(B(1,2))
       IF (AIMAG(B(1,2)) .EQ. D_ZERO) THEN
          R1 = B(1,1) / REAL(B(1,2))
       ELSE ! B(1,2) complex
          R1 = ZZDIV(B(1,1), B(1,2))
       END IF
       IF (ABS(B(2,2)) .GE. ABS(B(2,1))) THEN
          IF (AIMAG(B(2,2)) .EQ. D_ZERO) THEN
             R2 = B(2,1) / REAL(B(2,2))
          ELSE ! B(2,2) complex
             R2 = ZZDIV(B(2,1), B(2,2))
          END IF
          DO J = 1, N
             XX = ZZFMA(R1, X(I), Y(I)) !R1 * X(I) + Y(I)
             YY = ZZFMA(R2, X(I), Y(I)) !R2 * X(I) + Y(I)
             X(I) = ZZMUL(XX, B(1,2))
             Y(I) = ZZMUL(YY, B(2,2))
             I = I + LDA
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(2,1))
          IF (AIMAG(B(2,1)) .EQ. D_ZERO) THEN
             R2 = B(2,2) / REAL(B(2,1))
          ELSE ! B(2,1) complex
             R2 = ZZDIV(B(2,2), B(2,1))
          END IF
          DO J = 1, N
             XX = ZZFMA(R1, X(I), Y(I)) !R1 * X(I) + Y(I)
             YY = ZZFMA(R2, Y(I), X(I)) !X(I) + R2 * Y(I)
             X(I) = ZZMUL(XX, B(1,2))
             Y(I) = ZZMUL(YY, B(2,1))
             I = I + LDA
          END DO
       END IF
    END IF
  END SUBROUTINE BA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE AC(B, M, X, Y)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(IN) :: B(2,2)
    INTEGER, INTENT(IN) :: M
    COMPLEX(KIND=DWP), INTENT(INOUT) :: X(M), Y(M)

    COMPLEX(KIND=DWP) :: XX, YY
    INTEGER :: I

    ! DO I = 1, M
    !    XX = X(I) + Y(I) * B(2,1)
    !    YY = X(I) * B(1,2) + Y(I)
    !    X(I) = XX * REAL(B(1,1))
    !    Y(I) = YY * REAL(B(2,2))
    ! END DO

    DO I = 1, M
       XX = ZFMA(Y(I), B(2,1), X(I))
       YY = ZFMA(X(I), B(1,2), Y(I))
       X(I) = XX * REAL(B(1,1))
       Y(I) = YY * REAL(B(2,2))
    END DO
  END SUBROUTINE AC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE AB(B, M, X, Y)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(IN) :: B(2,2)
    INTEGER, INTENT(IN) :: M
    COMPLEX(KIND=DWP), INTENT(INOUT) :: X(M), Y(M)

    COMPLEX(KIND=DWP) :: R1, R2, XX, YY
    INTEGER :: I

    ! DO I = 1, M
    !    XX = X(I) * B(1,1) + Y(I) * B(2,1)
    !    YY = X(I) * B(1,2) + Y(I) * B(2,2)
    !    X(I) = XX
    !    Y(I) = YY
    ! END DO

    IF (ABS(B(1,1)) .GE. ABS(B(2,1))) THEN
       IF (AIMAG(B(1,1)) .EQ. D_ZERO) THEN
          R1 = B(2,1) / REAL(B(1,1))
       ELSE ! B(1,1) complex
          R1 = ZZDIV(B(2,1), B(1,1))
       END IF
       IF (ABS(B(2,2)) .GE. ABS(B(1,2))) THEN
          IF (AIMAG(B(2,2)) .EQ. D_ZERO) THEN
             R2 = B(1,2) / REAL(B(2,2))
          ELSE ! B(2,2) complex
             R2 = ZZDIV(B(1,2), B(2,2))
          END IF
          DO I = 1, M
             XX = ZZFMA(Y(I), R1, X(I)) !X(I) + Y(I) * R1
             YY = ZZFMA(X(I), R2, Y(I)) !X(I) * R2 + Y(I)
             X(I) = ZZMUL(XX, B(1,1))
             Y(I) = ZZMUL(YY, B(2,2))
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(1,2))
          IF (AIMAG(B(1,2)) .EQ. D_ZERO) THEN
             R2 = B(2,2) / REAL(B(1,2))
          ELSE ! B(1,2) complex
             R2 = ZZDIV(B(2,2), B(1,2))
          END IF
          DO I = 1, M
             XX = ZZFMA(Y(I), R1, X(I)) !X(I) + Y(I) * R1
             YY = ZZFMA(Y(I), R2, X(I)) !X(I) + Y(I) * R2
             X(I) = ZZMUL(XX, B(1,1))
             Y(I) = ZZMUL(YY, B(1,2))
          END DO
       END IF
    ELSE ! ABS(B(1,1)) .LT. ABS(B(2,1))
       IF (AIMAG(B(2,1)) .EQ. D_ZERO) THEN
          R1 = B(1,1) / REAL(B(2,1))
       ELSE ! B(2,1) complex
          R1 = ZZDIV(B(1,1), B(2,1))
       END IF
       IF (ABS(B(2,2)) .GE. ABS(B(1,2))) THEN
          IF (AIMAG(B(2,2)) .EQ. D_ZERO) THEN
             R2 = B(1,2) / REAL(B(2,2))
          ELSE ! B(2,2) complex
             R2 = ZZDIV(B(1,2), B(2,2))
          END IF
          DO I = 1, M
             XX = ZZFMA(X(I), R1, Y(I)) !X(I) * R1 + Y(I)
             YY = ZZFMA(X(I), R2, Y(I)) !X(I) * R2 + Y(I)
             X(I) = ZZMUL(XX, B(2,1))
             Y(I) = ZZMUL(YY, B(2,2))
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(1,2))
          IF (AIMAG(B(1,2)) .EQ. D_ZERO) THEN
             R2 = B(2,2) / REAL(B(1,2))
          ELSE ! B(1,2) complex
             R2 = ZZDIV(B(2,2), B(1,2))
          END IF
          DO I = 1, M
             XX = ZZFMA(X(I), R1, Y(I)) !X(I) * R1 + Y(I)
             YY = ZZFMA(Y(I), R2, X(I)) !X(I) + Y(I) * R2
             X(I) = ZZMUL(XX, B(2,1))
             Y(I) = ZZMUL(YY, B(1,2))
          END DO
       END IF
    END IF
  END SUBROUTINE AB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE REAL(KIND=DWP) FUNCTION ABODNZF2(B, M, X, Y, P, Q)
    ! skipping the diagonal elements X(P) and Y(Q),
    ! Z = SUM(|oldX(I)|^2-|newX(I)|^2 + |oldY(I)|^2-|newY(I)|^2)
    !   = ||oldX oldY||_F^2 - ||newX newY||_F^2
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: M, P, Q
    COMPLEX(KIND=DWP), INTENT(IN) :: B(2,2), X(M), Y(M)

    COMPLEX(KIND=DWP) :: R1, R2, XX, YY
    REAL(KIND=DWP) :: AXI, AXX, AYI, AYY
    INTEGER :: I

    IF ((M .LT. 2) .OR. (P .LE. 0) .OR. (P .GT. M) .OR. (Q .LE. P) .OR. (Q .GT. M)) THEN
       ABODNZF2 = D_MZERO
       RETURN
    END IF
    ABODNZF2 = D_ZERO
    ABODNZF2 = ABODNZF2 + REAL(Y(P)) * REAL(Y(P))
    ABODNZF2 = ABODNZF2 + AIMAG(Y(P)) * AIMAG(Y(P))
    ABODNZF2 = ABODNZF2 + REAL(X(Q)) * REAL(X(Q))
    ABODNZF2 = ABODNZF2 + AIMAG(X(Q)) * AIMAG(X(Q))

    IF (ABS(B(1,1)) .GE. ABS(B(2,1))) THEN
       IF (AIMAG(B(1,1)) .EQ. D_ZERO) THEN
          R1 = B(2,1) / REAL(B(1,1))
       ELSE ! B(1,1) complex
          R1 = ZZDIV(B(2,1), B(1,1))
       END IF
       IF (ABS(B(2,2)) .GE. ABS(B(1,2))) THEN
          IF (AIMAG(B(2,2)) .EQ. D_ZERO) THEN
             R2 = B(1,2) / REAL(B(2,2))
          ELSE ! B(2,2) complex
             R2 = ZZDIV(B(1,2), B(2,2))
          END IF
          DO I = 1, P-1
             XX = ZZMUL(ZZFMA(Y(I), R1, X(I)), B(1,1)) !(X(I) + Y(I) * R1) * B(1,1)
             YY = ZZMUL(ZZFMA(X(I), R2, Y(I)), B(2,2)) !(X(I) * R2 + Y(I)) * B(2,2)
             AXI = ABS(X(I))
             AXX = ABS(XX)
             AYI = ABS(Y(I))
             AYY = ABS(YY)
             ABODNZF2 = ABODNZF2 + (AXI - AXX) * (AXI + AXX)
             ABODNZF2 = ABODNZF2 + (AYI - AYY) * (AYI + AYY)
          END DO
          DO I = P+1, Q-1
             XX = ZZMUL(ZZFMA(Y(I), R1, X(I)), B(1,1)) !(X(I) + Y(I) * R1) * B(1,1)
             YY = ZZMUL(ZZFMA(X(I), R2, Y(I)), B(2,2)) !(X(I) * R2 + Y(I)) * B(2,2)
             AXI = ABS(X(I))
             AXX = ABS(XX)
             AYI = ABS(Y(I))
             AYY = ABS(YY)
             ABODNZF2 = ABODNZF2 + (AXI - AXX) * (AXI + AXX)
             ABODNZF2 = ABODNZF2 + (AYI - AYY) * (AYI + AYY)
          END DO
          DO I = Q+1, M
             XX = ZZMUL(ZZFMA(Y(I), R1, X(I)), B(1,1)) !(X(I) + Y(I) * R1) * B(1,1)
             YY = ZZMUL(ZZFMA(X(I), R2, Y(I)), B(2,2)) !(X(I) * R2 + Y(I)) * B(2,2)
             AXI = ABS(X(I))
             AXX = ABS(XX)
             AYI = ABS(Y(I))
             AYY = ABS(YY)
             ABODNZF2 = ABODNZF2 + (AXI - AXX) * (AXI + AXX)
             ABODNZF2 = ABODNZF2 + (AYI - AYY) * (AYI + AYY)
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(1,2))
          IF (AIMAG(B(1,2)) .EQ. D_ZERO) THEN
             R2 = B(2,2) / REAL(B(1,2))
          ELSE ! B(1,2) complex
             R2 = ZZDIV(B(2,2), B(1,2))
          END IF
          DO I = 1, P-1
             XX = ZZMUL(ZZFMA(Y(I), R1, X(I)), B(1,1)) !(X(I) + Y(I) * R1) * B(1,1)
             YY = ZZMUL(ZZFMA(Y(I), R2, X(I)), B(1,2)) !(X(I) + Y(I) * R2) * B(1,2)
             AXI = ABS(X(I))
             AXX = ABS(XX)
             AYI = ABS(Y(I))
             AYY = ABS(YY)
             ABODNZF2 = ABODNZF2 + (AXI - AXX) * (AXI + AXX)
             ABODNZF2 = ABODNZF2 + (AYI - AYY) * (AYI + AYY)
          END DO
          DO I = P+1, Q-1
             XX = ZZMUL(ZZFMA(Y(I), R1, X(I)), B(1,1)) !(X(I) + Y(I) * R1) * B(1,1)
             YY = ZZMUL(ZZFMA(Y(I), R2, X(I)), B(1,2)) !(X(I) + Y(I) * R2) * B(1,2)
             AXI = ABS(X(I))
             AXX = ABS(XX)
             AYI = ABS(Y(I))
             AYY = ABS(YY)
             ABODNZF2 = ABODNZF2 + (AXI - AXX) * (AXI + AXX)
             ABODNZF2 = ABODNZF2 + (AYI - AYY) * (AYI + AYY)
          END DO
          DO I = Q+1, M
             XX = ZZMUL(ZZFMA(Y(I), R1, X(I)), B(1,1)) !(X(I) + Y(I) * R1) * B(1,1)
             YY = ZZMUL(ZZFMA(Y(I), R2, X(I)), B(1,2)) !(X(I) + Y(I) * R2) * B(1,2)
             AXI = ABS(X(I))
             AXX = ABS(XX)
             AYI = ABS(Y(I))
             AYY = ABS(YY)
             ABODNZF2 = ABODNZF2 + (AXI - AXX) * (AXI + AXX)
             ABODNZF2 = ABODNZF2 + (AYI - AYY) * (AYI + AYY)
          END DO
       END IF
    ELSE ! ABS(B(1,1)) .LT. ABS(B(2,1))
       IF (AIMAG(B(2,1)) .EQ. D_ZERO) THEN
          R1 = B(1,1) / REAL(B(2,1))
       ELSE ! B(2,1) complex
          R1 = ZZDIV(B(1,1), B(2,1))
       END IF
       IF (ABS(B(2,2)) .GE. ABS(B(1,2))) THEN
          IF (AIMAG(B(2,2)) .EQ. D_ZERO) THEN
             R2 = B(1,2) / REAL(B(2,2))
          ELSE ! B(2,2) complex
             R2 = ZZDIV(B(1,2), B(2,2))
          END IF
          DO I = 1, P-1
             XX = ZZMUL(ZZFMA(X(I), R1, Y(I)), B(2,1)) !(X(I) * R1 + Y(I)) * B(2,1)
             YY = ZZMUL(ZZFMA(X(I), R2, Y(I)), B(2,2)) !(X(I) * R2 + Y(I)) * B(2,2)
             AXI = ABS(X(I))
             AXX = ABS(XX)
             AYI = ABS(Y(I))
             AYY = ABS(YY)
             ABODNZF2 = ABODNZF2 + (AXI - AXX) * (AXI + AXX)
             ABODNZF2 = ABODNZF2 + (AYI - AYY) * (AYI + AYY)
          END DO
          DO I = P+1, Q-1
             XX = ZZMUL(ZZFMA(X(I), R1, Y(I)), B(2,1)) !(X(I) * R1 + Y(I)) * B(2,1)
             YY = ZZMUL(ZZFMA(X(I), R2, Y(I)), B(2,2)) !(X(I) * R2 + Y(I)) * B(2,2)
             AXI = ABS(X(I))
             AXX = ABS(XX)
             AYI = ABS(Y(I))
             AYY = ABS(YY)
             ABODNZF2 = ABODNZF2 + (AXI - AXX) * (AXI + AXX)
             ABODNZF2 = ABODNZF2 + (AYI - AYY) * (AYI + AYY)
          END DO
          DO I = Q+1, M
             XX = ZZMUL(ZZFMA(X(I), R1, Y(I)), B(2,1)) !(X(I) * R1 + Y(I)) * B(2,1)
             YY = ZZMUL(ZZFMA(X(I), R2, Y(I)), B(2,2)) !(X(I) * R2 + Y(I)) * B(2,2)
             AXI = ABS(X(I))
             AXX = ABS(XX)
             AYI = ABS(Y(I))
             AYY = ABS(YY)
             ABODNZF2 = ABODNZF2 + (AXI - AXX) * (AXI + AXX)
             ABODNZF2 = ABODNZF2 + (AYI - AYY) * (AYI + AYY)
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(1,2))
          IF (AIMAG(B(1,2)) .EQ. D_ZERO) THEN
             R2 = B(2,2) / REAL(B(1,2))
          ELSE ! B(1,2) complex
             R2 = ZZDIV(B(2,2), B(1,2))
          END IF
          DO I = 1, P-1
             XX = ZZMUL(ZZFMA(X(I), R1, Y(I)), B(2,1)) !(X(I) * R1 + Y(I)) * B(2,1)
             YY = ZZMUL(ZZFMA(Y(I), R2, X(I)), B(1,2)) !(X(I) + Y(I) * R2) * B(1,2)
             AXI = ABS(X(I))
             AXX = ABS(XX)
             AYI = ABS(Y(I))
             AYY = ABS(YY)
             ABODNZF2 = ABODNZF2 + (AXI - AXX) * (AXI + AXX)
             ABODNZF2 = ABODNZF2 + (AYI - AYY) * (AYI + AYY)
          END DO
          DO I = P+1, Q-1
             XX = ZZMUL(ZZFMA(X(I), R1, Y(I)), B(2,1)) !(X(I) * R1 + Y(I)) * B(2,1)
             YY = ZZMUL(ZZFMA(Y(I), R2, X(I)), B(1,2)) !(X(I) + Y(I) * R2) * B(1,2)
             AXI = ABS(X(I))
             AXX = ABS(XX)
             AYI = ABS(Y(I))
             AYY = ABS(YY)
             ABODNZF2 = ABODNZF2 + (AXI - AXX) * (AXI + AXX)
             ABODNZF2 = ABODNZF2 + (AYI - AYY) * (AYI + AYY)
          END DO
          DO I = Q+1, M
             XX = ZZMUL(ZZFMA(X(I), R1, Y(I)), B(2,1)) !(X(I) * R1 + Y(I)) * B(2,1)
             YY = ZZMUL(ZZFMA(Y(I), R2, X(I)), B(1,2)) !(X(I) + Y(I) * R2) * B(1,2)
             AXI = ABS(X(I))
             AXX = ABS(XX)
             AYI = ABS(Y(I))
             AYY = ABS(YY)
             ABODNZF2 = ABODNZF2 + (AXI - AXX) * (AXI + AXX)
             ABODNZF2 = ABODNZF2 + (AYI - AYY) * (AYI + AYY)
          END DO
       END IF
    END IF
  END FUNCTION ABODNZF2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE ZHSVD2D(A, U, INFO)
    ! A diagonal
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    COMPLEX(KIND=DWP) :: V
    REAL(KIND=DWP) :: W

    INFO = 0
    A(2,1) = Z_ZERO
    A(1,2) = Z_ZERO

    IF (AIMAG(A(1,1)) .EQ. D_ZERO) THEN
       ! A(1,1) real
       IF (SIGN(D_ONE, REAL(A(1,1))) .EQ. D_MONE) THEN
          ! A(1,1) negative
          U(1,1) = -U(1,1)
          U(1,2) = -U(1,2)
          A(1,1) = -A(1,1)
       END IF
    ELSE IF (REAL(A(1,1)) .EQ. D_ZERO) THEN
       ! A(1,1) imaginary .NE. 0
       IF (AIMAG(A(1,1)) .LT. D_ZERO) THEN
          ! A(1,1) = i * negative
          ! V = CMPLX(D_ZERO, D_ONE, DWP)
          U(1,1) = CMPLX(-AIMAG(U(1,1)), REAL(U(1,1)), DWP)
          U(1,2) = CMPLX(-AIMAG(U(1,2)), REAL(U(1,2)), DWP)
          A(1,1) = CMPLX(-AIMAG(A(1,1)), REAL(A(1,1)), DWP)
       ELSE
          ! A(1,1) = i * positive
          ! V = CMPLX(D_ZERO, D_MONE, DWP)
          U(1,1) = CMPLX(AIMAG(U(1,1)), -REAL(U(1,1)), DWP)
          U(1,2) = CMPLX(AIMAG(U(1,2)), -REAL(U(1,2)), DWP)
          A(1,1) = CMPLX(AIMAG(A(1,1)), -REAL(A(1,1)), DWP)
       END IF
    ELSE
       ! A(1,1) complex .NE. 0
       W = ABS(A(1,1))
       V = CONJG(A(1,1) / W)
       U(1,1) = ZZMUL(V, U(1,1))
       U(1,2) = ZZMUL(V, U(1,2))
       A(1,1) = W
    END IF

    IF (AIMAG(A(2,2)) .EQ. D_ZERO) THEN
       ! A(2,2) real
       IF (SIGN(D_ONE, REAL(A(2,2))) .EQ. D_MONE) THEN
          ! A(2,2) negative
          U(2,1) = -U(2,1)
          U(2,2) = -U(2,2)
          A(2,2) = -A(2,2)
       END IF
    ELSE IF (REAL(A(2,2)) .EQ. D_ZERO) THEN
       ! A(2,2) imaginary .NE. 0
       IF (AIMAG(A(2,2)) .LT. D_ZERO) THEN
          ! A(2,2) = i * negative
          ! V = CMPLX(D_ZERO, D_ONE, DWP)
          U(2,1) = CMPLX(-AIMAG(U(2,1)), REAL(U(2,1)), DWP)
          U(2,2) = CMPLX(-AIMAG(U(2,2)), REAL(U(2,2)), DWP)
          A(2,2) = CMPLX(-AIMAG(A(2,2)), REAL(A(2,2)), DWP)
       ELSE
          ! A(2,2) = i * positive
          ! V = CMPLX(D_ZERO, D_MONE, DWP)
          U(2,1) = CMPLX(AIMAG(U(2,1)), -REAL(U(2,1)), DWP)
          U(2,2) = CMPLX(AIMAG(U(2,2)), -REAL(U(2,2)), DWP)
          A(2,2) = CMPLX(AIMAG(A(2,2)), -REAL(A(2,2)), DWP)
       END IF
    ELSE
       ! A(2,2) complex .NE. 0
       W = ABS(A(2,2))
       V = CONJG(A(2,2) / W)
       U(2,1) = ZZMUL(V, U(2,1))
       U(2,2) = ZZMUL(V, U(2,2))
       A(2,2) = W
    END IF
  END SUBROUTINE ZHSVD2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE ZHSVD2U(H, A, U, Z, INFO)
    ! A upper triangular, not diagonal
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: H
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    COMPLEX(KIND=DWP) :: V(2,2), W(2,2), X_, Y_, Z_
    REAL(KIND=DWP) :: X, Y, T2, TU, TZ, CU, CZ

    INFO = 0

    V(1,1) = Z_ONE
    V(2,1) = Z_ZERO
    V(1,2) = Z_ZERO
    V(2,2) = Z_ONE

    W(1,1) = Z_ONE
    W(2,1) = Z_ZERO
    W(1,2) = Z_ZERO
    W(2,2) = Z_ONE

    IF (H) THEN
       ! TODO: transform
       CONTINUE
    ELSE ! trigonometric
       X_ = CONJG(A(1,2) / REAL(A(1,1)))
       X = ABS(X_)
       Y = REAL(A(2,2)) / REAL(A(1,1))
       T2 = D_ONE + (X - Y) * (X + Y)
       IF (X .LT. Y) THEN
          T2 = -(SCALE(X, 1) * Y) / T2
       ELSE ! .GE.
          T2 = -(SCALE(Y, 1) * X) / T2
       END IF
       TU = T2 / (D_ONE + SQRT(D_ONE + T2 * T2))
       CU = D_ONE / SQRT(D_ONE + TU * TU)
       Y_ = A(1,2) * REAL(A(2,2))
       IF (Y_ .EQ. Z_ZERO) THEN
          Y_ = TU * SIGN(D_ONE, REAL(Y_))
       ELSE ! .NE. ZERO
          Y_ = (CONJG(Y_) / ABS(Y_)) * TU * SIGN(D_ONE, REAL(Y_))
       END IF
       V(1,1) = CU
       V(2,1) = Y_
       V(1,2) = -CONJG(Y_)
       V(2,2) = CU
       Z_ = ZDFMA(Y, Y_, -X_)
       TZ = ABS(Z_)
       CZ = D_ONE / SQRT(D_ONE + TZ * TZ)
       W(1,1) = CZ
       W(2,1) = -Z_
       W(1,2) = CONJG(Z_)
       W(2,2) = CZ
    END IF

    CALL CA(V, 2, U(1,1), U(2,1), 2)
    CALL CA(V, 2, A(1,1), A(2,1), 2)

    CALL AC(W, 2, A(1,1), A(1,2))
    CALL AC(W, 2, Z(1,1), Z(1,2))

    CALL ZHSVD2D(A, U, INFO)
    INFO = 1
  END SUBROUTINE ZHSVD2U

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE ZHSVD2L(A, U, Z, INFO)
    ! A lower triangular, not diagonal, hyperbolic J
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    COMPLEX(KIND=DWP) :: V(2,2), W(2,2)

    INFO = 0

    V(1,1) = Z_ONE
    V(2,1) = Z_ZERO
    V(1,2) = Z_ZERO
    V(2,2) = Z_ONE

    W(1,1) = Z_ONE
    W(2,1) = Z_ZERO
    W(1,2) = Z_ZERO
    W(2,2) = Z_ONE

    ! TODO: transform

    CALL CA(V, 2, U(1,1), U(2,1), 2)
    CALL CA(V, 2, A(1,1), A(2,1), 2)

    CALL AC(W, 2, A(1,1), A(1,2))
    CALL AC(W, 2, Z(1,1), Z(1,2))

    CALL ZHSVD2D(A, U, INFO)
    INFO = 1
  END SUBROUTINE ZHSVD2L

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE ZHSVD2G(H, A, U, Z, INFO)
    ! A general
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: H
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    REAL(KIND=DWP) :: C, T
    COMPLEX(KIND=DWP) :: Q(2,2), R, S

    ! REAL(KIND=DWP), EXTERNAL :: DZNRM2
    ! EXTERNAL :: ZLARTG, ZROT, ZSWAP

    R = CMPLX(ABS(A(1,1)), ABS(A(2,1)), DWP)
    S = CMPLX(ABS(A(1,2)), ABS(A(2,2)), DWP)
    C = ABS(R) ! ||A_1||
    T = ABS(S) ! ||A_2||

    ! column pivoting
    ! IF (DZNRM2(2, A(1,1), 1) .LT. DZNRM2(2, A(1,2), 1)) THEN
    IF (C .LT. T) THEN
       ! swap the columns of A
       ! CALL ZSWAP(2, A(1,1), 1, A(1,2), 1)
       R = A(1,1)
       S = A(2,1)
       A(1,1) = A(1,2)
       A(2,1) = A(2,2)
       A(1,2) = R
       A(2,2) = S
       ! swap the columns of Z
       ! IF (.NOT. H) CALL ZSWAP(2, Z(1,1), 1, Z(1,2), 1)
       IF (.NOT. H) THEN
          R = Z(1,1)
          S = Z(2,1)
          Z(1,1) = Z(1,2)
          Z(2,1) = Z(2,2)
          Z(1,2) = R
          Z(2,2) = S
       END IF
       ! record the swap
       INFO = 1
    ELSE ! no pivoting
       INFO = 0
    END IF

    ! row sorting
    IF (ABS(A(1,1)) .LT. ABS(A(2,1))) THEN
       ! swap the rows of U
       ! CALL ZSWAP(2, U(1,1), 2, U(2,1), 2)
       R = U(1,1)
       S = U(1,2)
       U(1,1) = U(2,1)
       U(1,2) = U(2,2)
       U(2,1) = R
       U(2,2) = S
       ! swap the rows of A
       ! CALL ZSWAP(2, A(1,1), 2, A(2,1), 2)
       R = A(1,1)
       S = A(1,2)
       A(1,1) = A(2,1)
       A(1,2) = A(2,2)
       A(2,1) = R
       A(2,2) = S
    END IF

    ! QR factorization of A
    ! CALL ZLARTG(A(1,1), A(2,1), C, S, R)
    ! IF (.NOT. (ABS(R) .LE. HUGE(R))) THEN
    !    INFO = -7
    !    RETURN
    ! END IF
    ! A(1,1) = R
    ! A(2,1) = Z_ZERO
    ! CALL ZROT(1, A(1,2), 2, A(2,2), 2, C, S)
    ! premultiply U by Q^H
    ! CALL ZROT(2, U(1,1), 2, U(2,1), 2, C, S)

    ! |A(1,1)| >= |A(2,1)| cannot be 0; otherwise,
    ! ||A_1||=0 >= ||A_2|| >= 0, so A=0
    ! specifically, A is diagonal

    ! C * |   1 e**i\alpha T | = Q
    !     | -e**-i\alpha T 1 |
    IF (AIMAG(A(1,1)) .EQ. D_ZERO) THEN
       S = A(2,1) / REAL(A(1,1))
    ELSE ! A(1,1) complex
       S = ZZDIV(A(2,1), A(1,1))
    END IF
    T = ABS(S)
    ! e**i\alpha = conjg(S)/T
    ! S is e**i\alpha * T here
    S = CONJG(S)
    C = D_ONE / SQRT(D_ONE + T * T)
    Q(1,1) =  C
    Q(2,1) = -CONJG(S)
    Q(1,2) =  S
    Q(2,2) =  C
    CALL CA(Q, 2, A(1,1), A(2,1), 2)
    ! should never happen
    IF (.NOT. (ABS(A(1,1)) .LE. HUGE(D_ZERO))) THEN
       INFO = -7
       RETURN
    END IF
    CALL CA(Q, 2, U(1,1), U(2,1), 2)
    A(2,1) = Z_ZERO

    ! make diag(A) real and non-negative
    ! the first row of U and A
    IF (AIMAG(A(1,1)) .EQ. D_ZERO) THEN
       ! A(1,1) real
       IF (SIGN(D_ONE, REAL(A(1,1))) .EQ. D_MONE) THEN
          ! A(1,1) negative
          U(1,1) = -U(1,1)
          U(1,2) = -U(1,2)
          A(1,1) = -A(1,1)
          A(1,2) = -A(1,2)
       END IF
    ELSE IF (REAL(A(1,1)) .EQ. D_ZERO) THEN
       ! A(1,1) imaginary .NE. 0
       IF (AIMAG(A(1,1)) .LT. D_ZERO) THEN
          ! A(1,1) = i * negative
          U(1,1) = CMPLX(-AIMAG(U(1,1)), REAL(U(1,1)), DWP)
          U(1,2) = CMPLX(-AIMAG(U(1,2)), REAL(U(1,2)), DWP)
          A(1,1) = CMPLX(-AIMAG(A(1,1)), REAL(A(1,1)), DWP)
          A(1,2) = CMPLX(-AIMAG(A(1,2)), REAL(A(1,2)), DWP)
       ELSE
          ! A(1,1) = i * positive
          U(1,1) = CMPLX(AIMAG(U(1,1)), -REAL(U(1,1)), DWP)
          U(1,2) = CMPLX(AIMAG(U(1,2)), -REAL(U(1,2)), DWP)
          A(1,1) = CMPLX(AIMAG(A(1,1)), -REAL(A(1,1)), DWP)
          A(1,2) = CMPLX(AIMAG(A(1,2)), -REAL(A(1,2)), DWP)
       END IF
    ELSE
       ! A(1,1) complex .NE. 0
       R = A(1,1)
       A(1,1) = ABS(A(1,1))
       R = CONJG(R / REAL(A(1,1)))
       U(1,1) = ZZMUL(R, U(1,1))
       U(1,2) = ZZMUL(R, U(1,2))
       A(1,2) = ZZMUL(R, A(1,2))
    END IF
    ! the second row of U and A
    IF (AIMAG(A(2,2)) .EQ. D_ZERO) THEN
       ! A(2,2) real
       IF (SIGN(D_ONE, REAL(A(2,2))) .EQ. D_MONE) THEN
          ! A(2,2) negative
          U(2,1) = -U(2,1)
          U(2,2) = -U(2,2)
          A(2,2) = -A(2,2)
       END IF
    ELSE IF (REAL(A(2,2)) .EQ. D_ZERO) THEN
       ! A(2,2) imaginary .NE. 0
       IF (AIMAG(A(2,2)) .LT. D_ZERO) THEN
          ! A(2,2) = i * negative
          U(2,1) = CMPLX(-AIMAG(U(2,1)), REAL(U(2,1)), DWP)
          U(2,2) = CMPLX(-AIMAG(U(2,2)), REAL(U(2,2)), DWP)
          A(2,2) = CMPLX(-AIMAG(A(2,2)), REAL(A(2,2)), DWP)
       ELSE
          ! A(2,2) = i * positive
          U(2,1) = CMPLX(AIMAG(U(2,1)), -REAL(U(2,1)), DWP)
          U(2,2) = CMPLX(AIMAG(U(2,2)), -REAL(U(2,2)), DWP)
          A(2,2) = CMPLX(AIMAG(A(2,2)), -REAL(A(2,2)), DWP)
       END IF
    ELSE
       ! A(2,2) complex .NE. 0
       R = A(2,2)
       A(2,2) = ABS(A(2,2))
       R = CONJG(R / REAL(A(2,2)))
       U(2,1) = ZZMUL(R, U(2,1))
       U(2,2) = ZZMUL(R, U(2,2))
    END IF

    IF (A(1,2) .EQ. Z_ZERO) THEN
       ! A diagonal
       IF (H) THEN
          IF (INFO .EQ. 1) THEN
             ! swap the rows of U
             ! CALL ZSWAP(2, U(1,1), 2, U(2,1), 2)
             R = U(1,1)
             S = U(1,2)
             U(1,1) = U(2,1)
             U(1,2) = U(2,2)
             U(2,1) = R
             U(2,2) = S
             ! swap the diagonal elements of A
             ! CALL ZSWAP(1, A(1,1), 1, A(2,2), 1)
             R = A(1,1)
             A(1,1) = A(2,2)
             A(2,2) = R
          ELSE ! mark the transform
             INFO = 1
          END IF
       ELSE ! mark the transform
          INFO = 1
       END IF
    ELSE IF (H .AND. (INFO .EQ. 1)) THEN
       ! swap the columns of A
       ! CALL ZSWAP(2, A(1,1), 1, A(1,2), 1)
       R = A(1,1)
       S = A(2,1)
       A(1,1) = A(1,2)
       A(2,1) = A(2,2)
       A(1,2) = R
       A(2,2) = S
       ! A upper antitriangular, not antidiagonal (X .NE. 0)
       !     | X R | <- R .NE. 0 the largest
       ! A = | x 0 |    element by magnitude
       ! swap the rows of U
       ! CALL ZSWAP(2, U(1,1), 2, U(2,1), 2)
       R = U(1,1)
       S = U(1,2)
       U(1,1) = U(2,1)
       U(1,2) = U(2,2)
       U(2,1) = R
       U(2,2) = S
       ! swap the rows of A
       ! CALL ZSWAP(2, A(1,1), 2, A(2,1), 2)
       R = A(1,1)
       S = A(1,2)
       A(1,1) = A(2,1)
       A(1,2) = A(2,2)
       A(2,1) = R
       A(2,2) = S
       ! A lower triangular, not diagonal (X .NE. 0)
       ! A = | x 0 |    R .NE. 0 the largest
       !     | X R | <- element by magnitude
       CALL ZHSVD2L(A, U, Z, INFO)
    ELSE
       ! A upper triangular, not diagonal (X .NE. 0)
       !     | R X | <- R .NE. 0 the largest
       ! A = | 0 x |    element by magnitude
       CALL ZHSVD2U(H, A, U, Z, INFO)
    END IF
  END SUBROUTINE ZHSVD2G

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE ZHSVD2S(H, A, U, Z, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: H
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    COMPLEX(KIND=DWP) :: R, S
    ! EXTERNAL :: ZSWAP

    ! assume A real, non-negative, diagonal
    ! permute A, U, Z in trigonometric case
    IF (((H .EQ. 2) .AND. (REAL(A(1,1)) .LT. REAL(A(2,2)))) .OR. ((H .EQ. -2) .AND. (REAL(A(1,1)) .GT. REAL(A(2,2))))) THEN
       ! swap the rows of U
       ! CALL ZSWAP(2, U(1,1), 2, U(2,1), 2)
       R = U(1,1)
       S = U(1,2)
       U(1,1) = U(2,1)
       U(1,2) = U(2,2)
       U(2,1) = R
       U(2,2) = S
       ! swap the diagonal elements of A
       ! CALL ZSWAP(1, A(1,1), 1, A(2,2), 1)
       R = A(1,1)
       A(1,1) = A(2,2)
       A(2,2) = R
       ! swap the columns of Z
       ! CALL ZSWAP(2, Z(1,1), 1, Z(1,2), 1)
       R = Z(1,1)
       S = Z(2,1)
       Z(1,1) = Z(1,2)
       Z(2,1) = Z(2,2)
       Z(1,2) = R
       Z(2,2) = S
    END IF

    ! check if U is identity and record in INFO if it is not
    IF ((U(1,1) .NE. Z_ONE) .OR. (U(2,1) .NE. Z_ZERO) .OR. (U(1,2) .NE. Z_ZERO) .OR. (U(2,2) .NE. Z_ONE)) INFO = INFO + 2
    ! check if Z is identity and record in INFO if it is not
    IF ((Z(1,1) .NE. Z_ONE) .OR. (Z(2,1) .NE. Z_ZERO) .OR. (Z(1,2) .NE. Z_ZERO) .OR. (Z(2,2) .NE. Z_ONE)) INFO = INFO + 4
  END SUBROUTINE ZHSVD2S

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE DSCALEW(W, S)
    IMPLICIT NONE

    REAL(KIND=DWP), PARAMETER :: TOOBIG = SCALE(HUGE(TOOBIG), -1)
    REAL(KIND=DWP), PARAMETER :: TOOSMALL = TINY(TOOSMALL)

    REAL(KIND=DWP), INTENT(INOUT) :: W(2,2)
    INTEGER, INTENT(OUT) :: S

    REAL(KIND=DWP) :: AA, AX
    INTEGER :: DS, US, I, J

    DS = 0
    US = 0

    DO J = 1, 2
       DO I = 1, 2
          AA = W(I,J)
          IF (.NOT. (AA .LE. HUGE(AA))) THEN
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
             W(I,J) = SCALE(W(I,J), DS)
          END DO
       END DO

       S = DS
    ELSE ! might perform upscaling
       AX = D_ZERO

       DO J = 1, 2
          DO I = 1, 2
             AA = W(I,J)
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
                W(I,J) = SCALE(W(I,J), US)
             END DO
          END DO
       END IF

       S = US
    END IF
  END SUBROUTINE DSCALEW

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE ZHSVD2(A, J, U, Z, INFO)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(2,2)
    INTEGER, INTENT(IN) :: J(2)
    COMPLEX(KIND=DWP), INTENT(OUT) :: U(2,2), Z(2,2)
    INTEGER, INTENT(OUT) :: INFO

    REAL(KIND=DWP) :: W(2,2)
    INTEGER :: S, I, K

    W(1,1) = ABS(A(1,1))
    W(2,1) = ABS(A(2,1))
    W(1,2) = ABS(A(1,2))
    W(2,2) = ABS(A(2,2))

    S = 0
    DO K = 1, 2
       DO I = 1, 2
          IF (.NOT. (W(I,K) .LE. HUGE(D_ZERO))) THEN
             IF ((ABS(REAL(A(I,K))) .LE. HUGE(D_ZERO)) .AND. (ABS(AIMAG(A(I,K))) .LE. HUGE(D_ZERO))) THEN
                S = -1
                GOTO 1
             END IF
             INFO = -((K - 1) * 2 + I)
             RETURN
          END IF
       END DO
    END DO

    ! prescale A
1   IF (S .NE. 0) THEN
       A(1,1) = CMPLX(SCALE(REAL(A(1,1)), S), SCALE(AIMAG(A(1,1)), S), DWP)
       A(2,1) = CMPLX(SCALE(REAL(A(2,1)), S), SCALE(AIMAG(A(2,1)), S), DWP)
       A(1,2) = CMPLX(SCALE(REAL(A(1,2)), S), SCALE(AIMAG(A(1,2)), S), DWP)
       A(2,2) = CMPLX(SCALE(REAL(A(2,2)), S), SCALE(AIMAG(A(2,2)), S), DWP)

       W(1,1) = ABS(A(1,1))
       W(2,1) = ABS(A(2,1))
       W(1,2) = ABS(A(1,2))
       W(2,2) = ABS(A(2,2))
    END IF
    K = S
    S = 0

    ! scale W as A would be scaled in the real case
    CALL DSCALEW(W, S)
    IF (S .LT. -1) THEN
       ! W has NaNs and/or infinities
       INFO = S + 1
       RETURN
    END IF

    SELECT CASE (J(1))
    CASE (-1,1)
       CONTINUE
    CASE DEFAULT
       INFO = -5
       RETURN
    END SELECT
    SELECT CASE (J(2))
    CASE (-1,1)
       CONTINUE
    CASE DEFAULT
       INFO = -6
       RETURN
    END SELECT
    INFO = 0

    ! scale A
    IF (S .NE. 0) THEN
       A(1,1) = CMPLX(SCALE(REAL(A(1,1)), S), SCALE(AIMAG(A(1,1)), S), DWP)
       A(2,1) = CMPLX(SCALE(REAL(A(2,1)), S), SCALE(AIMAG(A(2,1)), S), DWP)
       A(1,2) = CMPLX(SCALE(REAL(A(1,2)), S), SCALE(AIMAG(A(1,2)), S), DWP)
       A(2,2) = CMPLX(SCALE(REAL(A(2,2)), S), SCALE(AIMAG(A(2,2)), S), DWP)
       S = S + K
    ELSE ! S .EQ. 0
       S = K
    END IF

    ! U = I
    U(1,1) = Z_ONE
    U(2,1) = Z_ZERO
    U(1,2) = Z_ZERO
    U(2,2) = Z_ONE

    ! Z = I
    Z(1,1) = Z_ONE
    Z(2,1) = Z_ZERO
    Z(1,2) = Z_ZERO
    Z(2,2) = Z_ONE

    IF ((A(2,1) .EQ. Z_ZERO) .AND. (A(1,2) .EQ. Z_ZERO)) THEN
       ! A diagonal
       CALL ZHSVD2D(A, U, INFO)
    ELSE
       ! A general
       CALL ZHSVD2G((J(1) .NE. J(2)), A, U, Z, INFO)
    END IF
    IF (INFO .LT. 0) RETURN

    ! scale back if necessary
    IF (S .NE. 0) THEN
       I = -S
       A(1,1) = SCALE(REAL(A(1,1)), I)
       A(2,2) = SCALE(REAL(A(2,2)), I)
    END IF

    CALL ZHSVD2S((J(1) + J(2)), A, U, Z, INFO)
  END SUBROUTINE ZHSVD2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE ZTRANSF
