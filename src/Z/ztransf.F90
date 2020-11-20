MODULE ztransf
  USE ktransf
#ifndef ABSZ
#define ABSZ ABS
#endif
#ifndef FMAD
#define FMAD(a,b,c) ((a)*(b)+(c))
#endif
#ifndef MIND
#define MIND MIN
#endif
  IMPLICIT NONE

CONTAINS

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
    C = SIGN(MIND((ABS(R) / A), D_ONE), R)
    S = I / MAX(A, MINF)
    E = CMPLX(C, S, DWP)
  END SUBROUTINE ZPOLAR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ELEMENTAL COMPLEX(KIND=DWP) FUNCTION ZDFMA(A, B, C)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: A
    COMPLEX(KIND=DWP), INTENT(IN) :: B, C

    REAL(KIND=DWP) :: R, I

    R = FMAD(A, REAL(B), REAL(C))
    I = FMAD(A, AIMAG(B), AIMAG(C))

    ZDFMA = CMPLX(R, I, DWP)
  END FUNCTION ZDFMA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ELEMENTAL COMPLEX(KIND=DWP) FUNCTION ZJFMA(A, B, C)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: A ! imaginary
    COMPLEX(KIND=DWP), INTENT(IN) :: B, C

    REAL(KIND=DWP) :: R, I

    R = FMAD(-A, AIMAG(B), REAL(C))
    I = FMAD( A, REAL(B), AIMAG(C))

    ZJFMA = CMPLX(R, I, DWP)
  END FUNCTION ZJFMA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! after cuCfma() from CUDA's cuComplex.h

  ELEMENTAL COMPLEX(KIND=DWP) FUNCTION ZZFMA(A, B, C)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(IN) :: A, B, C

    REAL(KIND=DWP) :: F, R, I

    F = FMAD(-AIMAG(A), AIMAG(B), REAL(C))
    R = FMAD(REAL(A), REAL(B), F)
    F = FMAD(AIMAG(A), REAL(B), AIMAG(C))
    I = FMAD(REAL(A), AIMAG(B), F)

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
    R = FMAD(REAL(A), REAL(B), -F)
    F = AIMAG(A) * REAL(B)
    I = FMAD(REAL(A), AIMAG(B), F)

    ZZMUL = CMPLX(R, I, DWP)
  END FUNCTION ZZMUL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ELEMENTAL COMPLEX(KIND=DWP) FUNCTION ZDDIV(A, B)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(IN) :: A
    REAL(KIND=DWP), INTENT(IN) :: B

    ZDDIV = CMPLX(REAL(A) / B, AIMAG(A) / B, DWP)
  END FUNCTION ZDDIV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! assume that B is not real and |B| .NE. 0
  ELEMENTAL COMPLEX(KIND=DWP) FUNCTION ZZDIV(A, B)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(IN) :: A, B
    REAL(KIND=DWP) :: F

    F = ABSZ(B)
    ZZDIV = ZDDIV(ZZMUL(CONJG(ZDDIV(B, F)), A), F)
  END FUNCTION ZZDIV
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE C1A(B, A)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: B(2,2)
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(2)

    COMPLEX(KIND=DWP) :: C(2)

    C(1) = ZDDIV(ZDFMA(B(1,2), A(2), A(1)), B(1,1))
    C(2) = ZDDIV(ZDFMA(B(2,1), A(1), A(2)), B(2,2))
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
       C(1,J) = ZDDIV(ZDFMA(B(1,2), A(2,J), A(1,J)), B(1,1))
       C(2,J) = ZDDIV(ZDFMA(B(2,1), A(1,J), A(2,J)), B(2,2))
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
       C(I,1) = ZDDIV(ZDFMA(B(2,1), A(I,2), A(I,1)), B(1,1))
       C(I,2) = ZDDIV(ZDFMA(B(1,2), A(I,1), A(I,2)), B(2,2))
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

    IF (ABSZ(B(1,1)) .GE. ABSZ(B(1,2))) THEN
       IF (AIMAG(B(1,1)) .EQ. D_ZERO) THEN
          R1 = ZDDIV(B(1,2), REAL(B(1,1)))
       ELSE ! B(1,1) complex
          R1 = ZZDIV(B(1,2), B(1,1))
       END IF
       IF (ABSZ(B(2,2)) .GE. ABSZ(B(2,1))) THEN
          IF (AIMAG(B(2,2)) .EQ. D_ZERO) THEN
             R2 = ZDDIV(B(2,1), REAL(B(2,2)))
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
             R2 = ZDDIV(B(2,2), REAL(B(2,1)))
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
          R1 = ZDDIV(B(1,1), REAL(B(1,2)))
       ELSE ! B(1,2) complex
          R1 = ZZDIV(B(1,1), B(1,2))
       END IF
       IF (ABSZ(B(2,2)) .GE. ABSZ(B(2,1))) THEN
          IF (AIMAG(B(2,2)) .EQ. D_ZERO) THEN
             R2 = ZDDIV(B(2,1), REAL(B(2,2)))
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
             R2 = ZDDIV(B(2,2), REAL(B(2,1)))
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

    IF (ABSZ(B(1,1)) .GE. ABSZ(B(2,1))) THEN
       IF (AIMAG(B(1,1)) .EQ. D_ZERO) THEN
          R1 = ZDDIV(B(2,1), REAL(B(1,1)))
       ELSE ! B(1,1) complex
          R1 = ZZDIV(B(2,1), B(1,1))
       END IF
       IF (ABSZ(B(2,2)) .GE. ABSZ(B(1,2))) THEN
          IF (AIMAG(B(2,2)) .EQ. D_ZERO) THEN
             R2 = ZDDIV(B(1,2), REAL(B(2,2)))
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
             R2 = ZDDIV(B(2,2), REAL(B(1,2)))
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
          R1 = ZDDIV(B(1,1), REAL(B(2,1)))
       ELSE ! B(2,1) complex
          R1 = ZZDIV(B(1,1), B(2,1))
       END IF
       IF (ABSZ(B(2,2)) .GE. ABSZ(B(1,2))) THEN
          IF (AIMAG(B(2,2)) .EQ. D_ZERO) THEN
             R2 = ZDDIV(B(1,2), REAL(B(2,2)))
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
             R2 = ZDDIV(B(2,2), REAL(B(1,2)))
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

    ABODNZF2 = DASUM4(REAL(Y(P)), AIMAG(Y(P)), REAL(X(Q)), AIMAG(X(Q)))

    IF (ABSZ(B(1,1)) .GE. ABSZ(B(2,1))) THEN
       IF (AIMAG(B(1,1)) .EQ. D_ZERO) THEN
          R1 = ZDDIV(B(2,1), REAL(B(1,1)))
       ELSE ! B(1,1) complex
          R1 = ZZDIV(B(2,1), B(1,1))
       END IF
       IF (ABSZ(B(2,2)) .GE. ABSZ(B(1,2))) THEN
          IF (AIMAG(B(2,2)) .EQ. D_ZERO) THEN
             R2 = ZDDIV(B(1,2), REAL(B(2,2)))
          ELSE ! B(2,2) complex
             R2 = ZZDIV(B(1,2), B(2,2))
          END IF
          DO I = 1, P-1
             XX = ZZMUL(ZZFMA(Y(I), R1, X(I)), B(1,1)) !(X(I) + Y(I) * R1) * B(1,1)
             YY = ZZMUL(ZZFMA(X(I), R2, Y(I)), B(2,2)) !(X(I) * R2 + Y(I)) * B(2,2)
             AXI = ABSZ(X(I))
             AXX = ABSZ(XX)
             AYI = ABSZ(Y(I))
             AYY = ABSZ(YY)
             ABODNZF2 = FMAD((AXI - AXX), (AXI + AXX), ABODNZF2)
             ABODNZF2 = FMAD((AYI - AYY), (AYI + AYY), ABODNZF2)
          END DO
          DO I = P+1, Q-1
             XX = ZZMUL(ZZFMA(Y(I), R1, X(I)), B(1,1)) !(X(I) + Y(I) * R1) * B(1,1)
             YY = ZZMUL(ZZFMA(X(I), R2, Y(I)), B(2,2)) !(X(I) * R2 + Y(I)) * B(2,2)
             AXI = ABSZ(X(I))
             AXX = ABSZ(XX)
             AYI = ABSZ(Y(I))
             AYY = ABSZ(YY)
             ABODNZF2 = FMAD((AXI - AXX), (AXI + AXX), ABODNZF2)
             ABODNZF2 = FMAD((AYI - AYY), (AYI + AYY), ABODNZF2)
          END DO
          DO I = Q+1, M
             XX = ZZMUL(ZZFMA(Y(I), R1, X(I)), B(1,1)) !(X(I) + Y(I) * R1) * B(1,1)
             YY = ZZMUL(ZZFMA(X(I), R2, Y(I)), B(2,2)) !(X(I) * R2 + Y(I)) * B(2,2)
             AXI = ABSZ(X(I))
             AXX = ABSZ(XX)
             AYI = ABSZ(Y(I))
             AYY = ABSZ(YY)
             ABODNZF2 = FMAD((AXI - AXX), (AXI + AXX), ABODNZF2)
             ABODNZF2 = FMAD((AYI - AYY), (AYI + AYY), ABODNZF2)
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(1,2))
          IF (AIMAG(B(1,2)) .EQ. D_ZERO) THEN
             R2 = ZDDIV(B(2,2), REAL(B(1,2)))
          ELSE ! B(1,2) complex
             R2 = ZZDIV(B(2,2), B(1,2))
          END IF
          DO I = 1, P-1
             XX = ZZMUL(ZZFMA(Y(I), R1, X(I)), B(1,1)) !(X(I) + Y(I) * R1) * B(1,1)
             YY = ZZMUL(ZZFMA(Y(I), R2, X(I)), B(1,2)) !(X(I) + Y(I) * R2) * B(1,2)
             AXI = ABSZ(X(I))
             AXX = ABSZ(XX)
             AYI = ABSZ(Y(I))
             AYY = ABSZ(YY)
             ABODNZF2 = FMAD((AXI - AXX), (AXI + AXX), ABODNZF2)
             ABODNZF2 = FMAD((AYI - AYY), (AYI + AYY), ABODNZF2)
          END DO
          DO I = P+1, Q-1
             XX = ZZMUL(ZZFMA(Y(I), R1, X(I)), B(1,1)) !(X(I) + Y(I) * R1) * B(1,1)
             YY = ZZMUL(ZZFMA(Y(I), R2, X(I)), B(1,2)) !(X(I) + Y(I) * R2) * B(1,2)
             AXI = ABSZ(X(I))
             AXX = ABSZ(XX)
             AYI = ABSZ(Y(I))
             AYY = ABSZ(YY)
             ABODNZF2 = FMAD((AXI - AXX), (AXI + AXX), ABODNZF2)
             ABODNZF2 = FMAD((AYI - AYY), (AYI + AYY), ABODNZF2)
          END DO
          DO I = Q+1, M
             XX = ZZMUL(ZZFMA(Y(I), R1, X(I)), B(1,1)) !(X(I) + Y(I) * R1) * B(1,1)
             YY = ZZMUL(ZZFMA(Y(I), R2, X(I)), B(1,2)) !(X(I) + Y(I) * R2) * B(1,2)
             AXI = ABSZ(X(I))
             AXX = ABSZ(XX)
             AYI = ABSZ(Y(I))
             AYY = ABSZ(YY)
             ABODNZF2 = FMAD((AXI - AXX), (AXI + AXX), ABODNZF2)
             ABODNZF2 = FMAD((AYI - AYY), (AYI + AYY), ABODNZF2)
          END DO
       END IF
    ELSE ! ABS(B(1,1)) .LT. ABS(B(2,1))
       IF (AIMAG(B(2,1)) .EQ. D_ZERO) THEN
          R1 = ZDDIV(B(1,1), REAL(B(2,1)))
       ELSE ! B(2,1) complex
          R1 = ZZDIV(B(1,1), B(2,1))
       END IF
       IF (ABSZ(B(2,2)) .GE. ABSZ(B(1,2))) THEN
          IF (AIMAG(B(2,2)) .EQ. D_ZERO) THEN
             R2 = ZDDIV(B(1,2), REAL(B(2,2)))
          ELSE ! B(2,2) complex
             R2 = ZZDIV(B(1,2), B(2,2))
          END IF
          DO I = 1, P-1
             XX = ZZMUL(ZZFMA(X(I), R1, Y(I)), B(2,1)) !(X(I) * R1 + Y(I)) * B(2,1)
             YY = ZZMUL(ZZFMA(X(I), R2, Y(I)), B(2,2)) !(X(I) * R2 + Y(I)) * B(2,2)
             AXI = ABSZ(X(I))
             AXX = ABSZ(XX)
             AYI = ABSZ(Y(I))
             AYY = ABSZ(YY)
             ABODNZF2 = FMAD((AXI - AXX), (AXI + AXX), ABODNZF2)
             ABODNZF2 = FMAD((AYI - AYY), (AYI + AYY), ABODNZF2)
          END DO
          DO I = P+1, Q-1
             XX = ZZMUL(ZZFMA(X(I), R1, Y(I)), B(2,1)) !(X(I) * R1 + Y(I)) * B(2,1)
             YY = ZZMUL(ZZFMA(X(I), R2, Y(I)), B(2,2)) !(X(I) * R2 + Y(I)) * B(2,2)
             AXI = ABSZ(X(I))
             AXX = ABSZ(XX)
             AYI = ABSZ(Y(I))
             AYY = ABSZ(YY)
             ABODNZF2 = FMAD((AXI - AXX), (AXI + AXX), ABODNZF2)
             ABODNZF2 = FMAD((AYI - AYY), (AYI + AYY), ABODNZF2)
          END DO
          DO I = Q+1, M
             XX = ZZMUL(ZZFMA(X(I), R1, Y(I)), B(2,1)) !(X(I) * R1 + Y(I)) * B(2,1)
             YY = ZZMUL(ZZFMA(X(I), R2, Y(I)), B(2,2)) !(X(I) * R2 + Y(I)) * B(2,2)
             AXI = ABSZ(X(I))
             AXX = ABSZ(XX)
             AYI = ABSZ(Y(I))
             AYY = ABSZ(YY)
             ABODNZF2 = FMAD((AXI - AXX), (AXI + AXX), ABODNZF2)
             ABODNZF2 = FMAD((AYI - AYY), (AYI + AYY), ABODNZF2)
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(1,2))
          IF (AIMAG(B(1,2)) .EQ. D_ZERO) THEN
             R2 = ZDDIV(B(2,2), REAL(B(1,2)))
          ELSE ! B(1,2) complex
             R2 = ZZDIV(B(2,2), B(1,2))
          END IF
          DO I = 1, P-1
             XX = ZZMUL(ZZFMA(X(I), R1, Y(I)), B(2,1)) !(X(I) * R1 + Y(I)) * B(2,1)
             YY = ZZMUL(ZZFMA(Y(I), R2, X(I)), B(1,2)) !(X(I) + Y(I) * R2) * B(1,2)
             AXI = ABSZ(X(I))
             AXX = ABSZ(XX)
             AYI = ABSZ(Y(I))
             AYY = ABSZ(YY)
             ABODNZF2 = FMAD((AXI - AXX), (AXI + AXX), ABODNZF2)
             ABODNZF2 = FMAD((AYI - AYY), (AYI + AYY), ABODNZF2)
          END DO
          DO I = P+1, Q-1
             XX = ZZMUL(ZZFMA(X(I), R1, Y(I)), B(2,1)) !(X(I) * R1 + Y(I)) * B(2,1)
             YY = ZZMUL(ZZFMA(Y(I), R2, X(I)), B(1,2)) !(X(I) + Y(I) * R2) * B(1,2)
             AXI = ABSZ(X(I))
             AXX = ABSZ(XX)
             AYI = ABSZ(Y(I))
             AYY = ABSZ(YY)
             ABODNZF2 = FMAD((AXI - AXX), (AXI + AXX), ABODNZF2)
             ABODNZF2 = FMAD((AYI - AYY), (AYI + AYY), ABODNZF2)
          END DO
          DO I = Q+1, M
             XX = ZZMUL(ZZFMA(X(I), R1, Y(I)), B(2,1)) !(X(I) * R1 + Y(I)) * B(2,1)
             YY = ZZMUL(ZZFMA(Y(I), R2, X(I)), B(1,2)) !(X(I) + Y(I) * R2) * B(1,2)
             AXI = ABSZ(X(I))
             AXX = ABSZ(XX)
             AYI = ABSZ(Y(I))
             AYY = ABSZ(YY)
             ABODNZF2 = FMAD((AXI - AXX), (AXI + AXX), ABODNZF2)
             ABODNZF2 = FMAD((AYI - AYY), (AYI + AYY), ABODNZF2)
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

    D = .TRUE.
    ! C = ||A_1||
    IF (A(2,1) .EQ. Z_ZERO) THEN
       C = ABSZ(A(1,1))
    ELSE IF (A(1,1) .EQ. Z_ZERO) THEN
       C = ABSZ(A(2,1))
       D = .FALSE.
    ELSE ! full 1st column
       C = HYPOT(ABSZ(A(1,1)), ABSZ(A(2,1)))
       D = .FALSE.
    END IF
    ! T = ||A_2||
    IF (A(1,2) .EQ. Z_ZERO) THEN
       T = ABSZ(A(2,2))
    ELSE IF (A(2,2) .EQ. Z_ZERO) THEN
       T = ABSZ(A(1,2))
       D = .FALSE.
    ELSE ! full 2nd column
       T = HYPOT(ABSZ(A(1,2)), ABSZ(A(2,2)))
       D = .FALSE.
    END IF

    ! column pivoting
    IF (C .LT. T) THEN
       ! swap the columns of A
       R = A(1,1)
       S = A(2,1)
       A(1,1) = A(1,2)
       A(2,1) = A(2,2)
       A(1,2) = R
       A(2,2) = S
       ! swap the columns of Z
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
       R = U(1,1)
       S = U(1,2)
       U(1,1) = U(2,1)
       U(1,2) = U(2,2)
       U(2,1) = R
       U(2,2) = S
       ! swap the rows of A
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
       T = REAL(A(2,1)) / REAL(A(1,1))
       C = SQRT(FMAD(T, T, D_ONE)) ! D_ONE /
       Q(1,1) =  C
       Q(2,1) = -T
       Q(1,2) =  T
       Q(2,2) =  C

       ! premultiply U by Q^H
       CALL C2A(Q, U)
       CALL C1A(Q, A(1,2))
    END IF
    A(1,1) = CMPLX(W, D_ZERO, DWP)
    A(2,1) = Z_ZERO

    IF ((J(1) .NE. J(2)) .AND. P) THEN
       ! swap the columns of A
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
       R = U(1,1)
       S = U(1,2)
       U(1,1) = U(2,1)
       U(1,2) = U(2,2)
       U(2,1) = R
       U(2,2) = S
       ! swap the rows of A
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

    IF (((T .EQ. 2) .AND. (REAL(A(1,1)) .LT. REAL(A(2,2)))) .OR. ((T .EQ. -2) .AND. (REAL(A(2,2)) .LT. REAL(A(1,1))))) THEN
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
    ! IF (REAL(A(1,1)) .GT. HUGE(D_ZERO)) INFO = INFO + 8
    ! IF (REAL(A(2,2)) .GT. HUGE(D_ZERO)) INFO = INFO + 16
#endif
  END SUBROUTINE ZHSVD2S

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE DSCALEA(A, S)
    IMPLICIT NONE

    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(2,2)
    INTEGER, INTENT(OUT) :: S

    INTEGER :: I, J

    S = HUGE(S)
    DO J = 1, 2
       DO I = 1, 2
          IF (.NOT. (ABS(REAL(A(I,J))) .LE. HUGE(D_ZERO))) THEN
             ! -3 <= S <= -6
             S = -(J * 2 + I)
             RETURN
          END IF
          S = MIN(S, (EH - EXPONENT(REAL(A(I,J))) - 2))
          IF (.NOT. (ABS(AIMAG(A(I,J))) .LE. HUGE(D_ZERO))) THEN
             ! -3 <= S <= -6
             S = -(J * 2 + I)
             RETURN
          END IF
          S = MIN(S, (EH - EXPONENT(AIMAG(A(I,J))) - 2))
       END DO
    END DO

    IF (S .NE. 0) THEN
       DO J = 1, 2
          DO I = 1, 2
             A(I,J) = CMPLX(SCALE(REAL(A(I,J)),S), SCALE(AIMAG(A(I,J)),S), DWP)
          END DO
       END DO
    END IF
  END SUBROUTINE DSCALEA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
       INFO = -1
    ELSE IF (ABS(J(2)) .NE. 1) THEN ! error
       INFO = -2
    ELSE ! J OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN
#endif

    CALL DSCALEA(A, S)
    IF (S .LE. -3) THEN
       ! A has NaNs and/or infinities
       INFO = S
       RETURN
#ifdef NDEBUG
    ELSE
       INFO = 0
#endif
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
    IF (INFO .GT. 0) THEN
       W(1,1) = REAL(A(1,1))
       W(2,1) = REAL(A(2,1))
       W(1,2) = REAL(A(1,2))
       W(2,2) = REAL(A(2,2))
       CALL KHSVD2(W, J, B, C, S, INFO)
    ELSE IF (INFO .EQ. 0) THEN
       W(1,1) = SCALE(REAL(A(1,1)), -S)
       W(2,2) = SCALE(REAL(A(2,2)), -S)
    END IF
    IF (INFO .LT. 0) RETURN
    IF (INFO .GT. 0) THEN
       CALL C2A(B, U)
       CALL A2C(Z, C)
    END IF
    A(1,1) = CMPLX(W(1,1), D_ZERO, DWP)
    A(2,2) = CMPLX(W(2,2), D_ZERO, DWP)
    CALL ZHSVD2S(A, J, U, Z, INFO)
  END SUBROUTINE ZHSVD2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE ztransf
