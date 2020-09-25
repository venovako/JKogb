MODULE ZTRANSF
  USE KTRANSF
  IMPLICIT NONE

CONTAINS

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

  PURE SUBROUTINE ZPOLAR(Z, A, E)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(IN) :: Z
    REAL(KIND=DWP), INTENT(OUT) :: A
    COMPLEX(KIND=DWP), INTENT(OUT) :: E

    REAL(KIND=DWP) :: R, I, C, S

    R = REAL(Z)
    I = AIMAG(Z)
    A = HYPOT(R, I)
    C = SIGN(MIN(ABS(R) / A, D_ONE), R)
    S = I / MAX(A, MINF)
    E = CMPLX(C, S, DWP)
  END SUBROUTINE ZPOLAR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ELEMENTAL COMPLEX(KIND=DWP) FUNCTION ZDFMA(A, B, C)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: A
    COMPLEX(KIND=DWP), INTENT(IN) :: B, C

    REAL(KIND=DWP) :: R, I

    R = A * REAL(B) + REAL(C)
    I = A * AIMAG(B) + AIMAG(C)

    ZDFMA = CMPLX(R, I, DWP)
  END FUNCTION ZDFMA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ELEMENTAL COMPLEX(KIND=DWP) FUNCTION ZJFMA(A, B, C)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: A ! imaginary
    COMPLEX(KIND=DWP), INTENT(IN) :: B, C

    REAL(KIND=DWP) :: R, I

    R = REAL(C) - A * AIMAG(B)
    I = A * REAL(B) + AIMAG(C)

    ZJFMA = CMPLX(R, I, DWP)
  END FUNCTION ZJFMA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! after cuCfma() from CUDA's cuComplex.h

  ELEMENTAL COMPLEX(KIND=DWP) FUNCTION ZZFMA(A, B, C)
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

  ELEMENTAL COMPLEX(KIND=DWP) FUNCTION ZFMA(A, B, C)
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

  ELEMENTAL COMPLEX(KIND=DWP) FUNCTION ZZMUL(A, B)
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
  ELEMENTAL COMPLEX(KIND=DWP) FUNCTION ZZDIV(A, B)
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

  PURE SUBROUTINE C1A(B, A)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: B(2,2)
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(2)

    COMPLEX(KIND=DWP) :: C(2)

    C(1) = ZDFMA(B(1,2), A(2), A(1)) / B(1,1)
    C(2) = ZDFMA(B(2,1), A(1), A(2)) / B(2,2)
    A = C
  END SUBROUTINE C1A

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE C2A(B, A)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: B(2,2)
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(2,2)

    COMPLEX(KIND=DWP) :: C(2,2)
    INTEGER :: J

    DO J = 1, 2
       C(1,J) = ZDFMA(B(1,2), A(2,J), A(1,J)) / B(1,1)
       C(2,J) = ZDFMA(B(2,1), A(1,J), A(2,J)) / B(2,2)
    END DO
    A = C
  END SUBROUTINE C2A

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE A2C(A, B)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(2,2)
    REAL(KIND=DWP), INTENT(IN) :: B(2,2)

    COMPLEX(KIND=DWP) :: C(2,2)
    INTEGER :: I

    DO I = 1, 2
       C(I,1) = ZDFMA(B(2,1), A(I,2), A(I,1)) / B(1,1)
       C(I,2) = ZDFMA(B(1,2), A(I,1), A(I,2)) / B(2,2)
    END DO
    A = C
  END SUBROUTINE A2C

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

#ifndef NDEBUG
    IF ((M .LT. 2) .OR. (P .LE. 0) .OR. (P .GT. M) .OR. (Q .LE. P) .OR. (Q .GT. M)) THEN
       ABODNZF2 = D_MZERO
       RETURN
    END IF
#endif

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

  PURE SUBROUTINE ZHSVD2T(A, J, U, Z, INFO)
    ! A general
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(IN) :: J(2)
    INTEGER, INTENT(OUT) :: INFO

    REAL(KIND=DWP) :: Q(2,2), C, T, W
    COMPLEX(KIND=DWP) :: R, S
    LOGICAL :: P, D

    ! REAL(KIND=DWP), EXTERNAL :: DZNRM2
    ! EXTERNAL :: ZLARTG, ZROT, ZSWAP

    D = .TRUE.
    ! C = ||A_1||
    IF (A(2,1) .EQ. Z_ZERO) THEN
       C = ABS(A(1,1))
    ELSE IF (A(1,1) .EQ. Z_ZERO) THEN
       C = ABS(A(2,1))
       D = .FALSE.
    ELSE ! full 1st column
       C = HYPOT(ABS(A(1,1)), ABS(A(2,1)))
       D = .FALSE.
    END IF
    ! T = ||A_2||
    IF (A(1,2) .EQ. Z_ZERO) THEN
       T = ABS(A(2,2))
    ELSE IF (A(2,2) .EQ. Z_ZERO) THEN
       T = ABS(A(1,2))
       D = .FALSE.
    ELSE ! full 2nd column
       T = HYPOT(ABS(A(1,2)), ABS(A(2,2)))
       D = .FALSE.
    END IF

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
       IF (J(1) .EQ. J(2)) THEN
          R = Z(1,1)
          S = Z(2,1)
          Z(1,1) = Z(1,2)
          Z(2,1) = Z(2,2)
          Z(1,2) = R
          Z(2,2) = S
       END IF
       ! record the swap
       P = .TRUE.
       W = T
    ELSE ! no pivoting
       P = .FALSE.
       W = C
    END IF

    CALL ZPOLAR(A(1,1), C, R)
    ! make A(1,1) non-negative
    IF (R .NE. Z_ONE) THEN
       S = CONJG(R)
       U(1,1) = ZZMUL(U(1,1), S)
       U(1,2) = ZZMUL(U(1,2), S)
       A(1,1) = ZZMUL(A(1,1), S)
       A(1,2) = ZZMUL(A(1,2), S)
    END IF

    CALL ZPOLAR(A(2,1), T, S)
    ! make A(2,1) non-negative
    IF (S .NE. Z_ONE) THEN
       R = CONJG(S)
       U(2,1) = ZZMUL(U(2,1), R)
       U(2,2) = ZZMUL(U(2,2), R)
       A(2,1) = ZZMUL(A(2,1), R)
       A(2,2) = ZZMUL(A(2,2), R)
    END IF

    ! row sorting
    IF (REAL(A(1,1)) .LT. REAL(A(2,1))) THEN
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
    ! |A(1,1)| >= |A(2,1)| cannot be 0; otherwise,
    ! ||A_1||=0 >= ||A_2|| >= 0, so A=0
    ! specifically, A is diagonal

    IF (.NOT. D) THEN
       T = A(2,1) / A(1,1)
       C = SQRT(T * T + D_ONE) ! D_ONE /
       Q(1,1) =  C
       Q(2,1) = -T
       Q(1,2) =  T
       Q(2,2) =  C

       CALL C2A(Q, U)
       CALL C1A(Q, A(1,2))
    END IF
    A(1,1) = CMPLX(W, D_ZERO, DWP)
    A(2,1) = Z_ZERO

    IF ((J(1) .NE. J(2)) .AND. P) THEN
       ! swap the columns of A
       ! CALL ZSWAP(2, A(1,1), 1, A(1,2), 1)
       R = A(1,1)
       S = A(2,1)
       A(1,1) = A(1,2)
       A(2,1) = A(2,2)
       A(1,2) = R
       A(2,2) = S
       ! A upper antitriangular
       !     | x R | <- R is the largest
       ! A = | y 0 |    element by magnitude
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
       ! A lower triangular
       ! A = | y 0 |    R is the largest
       !     | x R | <- element by magnitude

       CALL ZPOLAR(A(2,1), T, S)
       ! make A(2,1) non-negative
       IF (S .NE. Z_ONE) THEN
          R = CONJG(S)
          Z(1,1) = ZZMUL(Z(1,1), R)
          Z(2,1) = ZZMUL(Z(2,1), R)
          A(1,1) = ZZMUL(A(1,1), R)
          A(2,1) = ZZMUL(A(2,1), R)
       END IF

       CALL ZPOLAR(A(1,1), C, R)
       ! make A(1,1) non-negative
       IF (R .NE. Z_ONE) THEN
          S = CONJG(R)
          U(1,1) = ZZMUL(U(1,1), S)
          U(1,2) = ZZMUL(U(1,2), S)
          A(1,1) = ZZMUL(A(1,1), S)
       END IF
    ELSE ! trigonometric or no column pivoting
       ! A upper triangular
       !     | R x | <- R is the largest
       ! A = | 0 y |    element by magnitude

       CALL ZPOLAR(A(1,2), C, R)
       ! make A(1,2) non-negative
       IF (R .NE. Z_ONE) THEN
          S = CONJG(R)
          Z(1,2) = ZZMUL(Z(1,2), S)
          Z(2,2) = ZZMUL(Z(2,2), S)
          A(1,2) = ZZMUL(A(1,2), S)
          A(2,2) = ZZMUL(A(2,2), S)
       END IF

       CALL ZPOLAR(A(2,2), T, S)
       ! make A(2,2) non-negative
       IF (S .NE. Z_ONE) THEN
          R = CONJG(S)
          U(2,1) = ZZMUL(U(2,1), R)
          U(2,2) = ZZMUL(U(2,2), R)
          A(2,2) = ZZMUL(A(2,2), R)
       END IF
    END IF

    IF (D) THEN
       INFO = 0
    ELSE ! .NOT. D
       INFO = 1
    END IF
  END SUBROUTINE ZHSVD2T

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE ZHSVD2S(A, J, U, Z, INFO)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(IN) :: J(2)
    INTEGER, INTENT(INOUT) :: INFO

    INTEGER :: T

    T = J(1) + J(2)
#ifndef NDEBUG
    IF (ABS(T) .NE. 2) INFO = -12
    IF (INFO .LT. 0) RETURN
#endif

    IF ((INFO .EQ. 1) .AND. &
         (((T .EQ. 2) .AND. (REAL(A(1,1)) .LT. REAL(A(2,2)))) .OR. ((T .EQ. -2) .AND. (REAL(A(2,2)) .LT. REAL(A(1,1)))))) THEN
       ! swap the rows of U
       A(2,1) = U(1,1)
       A(1,2) = U(1,2)
       U(1,1) = U(2,1)
       U(1,2) = U(2,2)
       U(2,1) = A(2,1)
       U(2,2) = A(1,2)
       ! swap the diagonal elements of A
       A(2,1) = A(1,1)
       A(1,1) = A(2,2)
       A(2,2) = A(2,1)
       ! swap the columns of Z
       A(2,1) = Z(1,1)
       A(1,2) = Z(2,1)
       Z(1,1) = Z(1,2)
       Z(2,1) = Z(2,2)
       Z(1,2) = A(2,1)
       Z(2,2) = A(1,2)
    END IF
    A(2,1) = Z_ZERO
    A(1,2) = Z_ZERO

    ! check if U is identity and record in INFO if it is not
    IF ((U(1,1) .NE. Z_ONE) .OR. (U(2,1) .NE. Z_ZERO) .OR. (U(1,2) .NE. Z_ZERO) .OR. (U(2,2) .NE. Z_ONE)) INFO = INFO + 2
    ! check if Z is identity and record in INFO if it is not
    IF ((Z(1,1) .NE. Z_ONE) .OR. (Z(2,1) .NE. Z_ZERO) .OR. (Z(1,2) .NE. Z_ZERO) .OR. (Z(2,2) .NE. Z_ONE)) INFO = INFO + 4
#ifndef NDEBUG
    ! check for overflow
    IF (REAL(A(1,1)) .GT. HUGE(D_ZERO)) INFO = INFO + 8
    IF (REAL(A(2,2)) .GT. HUGE(D_ZERO)) INFO = INFO + 16
#endif
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

! #ifndef USE_EXTENDED
  PURE SUBROUTINE ZHSVD2(A, J, U, Z, INFO)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(2,2)
    INTEGER, INTENT(IN) :: J(2)
    COMPLEX(KIND=DWP), INTENT(OUT) :: U(2,2), Z(2,2)
    INTEGER, INTENT(OUT) :: INFO

    REAL(KIND=DWP) :: W(2,2), B(2,2), C(2,2)
    INTEGER :: S, I, K

#ifndef NDEBUG
    IF (ABS(J(1)) .NE. 1) THEN ! error
       INFO = -5
    ELSE IF (ABS(J(2)) .NE. 1) THEN ! error
       INFO = -6
    ELSE ! J OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN
#endif

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
#ifdef NDEBUG
    ELSE
       INFO = 0
#endif
    END IF

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

    CALL ZHSVD2T(A, J, U, Z, INFO)
    W(1,1) = REAL(A(1,1))
    W(2,1) = REAL(A(2,1))
    W(1,2) = REAL(A(1,2))
    W(2,2) = REAL(A(2,2))
    CALL KHSVD2(W, J, B, C, INFO)
    IF (INFO .LT. 0) RETURN
    IF (INFO .EQ. 1) THEN
       CALL C2A(B, U)
       CALL A2C(Z, C)
    END IF
    IF (S .NE. 0) THEN
       A(1,1) = CMPLX(SCALE(W(1,1), -S), D_ZERO, DWP)
       A(2,2) = CMPLX(SCALE(W(2,2), -S), D_ZERO, DWP)
    ELSE ! no scaling
       A(1,1) = CMPLX(W(1,1), D_ZERO, DWP)
       A(2,2) = CMPLX(W(2,2), D_ZERO, DWP)
    END IF
    CALL ZHSVD2S(A, J, U, Z, INFO)
  END SUBROUTINE ZHSVD2
! #endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE ZTRANSF
