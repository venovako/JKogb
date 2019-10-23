MODULE ZTRANSF
  USE PARAMS
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

  PURE SUBROUTINE CA(B, N, X, Y, LDA)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(IN) :: B(2,2)
    INTEGER, INTENT(IN) :: N, LDA
    COMPLEX(KIND=DWP), INTENT(INOUT) :: X(*), Y(*)

    COMPLEX(KIND=DWP) :: XX, YY
    INTEGER :: I, J

    I = 1
    !DIR$ VECTOR ALWAYS
    DO J = 1, N
       XX = X(I) + B(1,2) * Y(I)
       YY = B(2,1) * X(I) + Y(I)
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
       R1 = B(1,2) / B(1,1)
       IF (ABS(B(2,2)) .GE. ABS(B(2,1))) THEN
          R2 = B(2,1) / B(2,2)
          !DIR$ VECTOR ALWAYS
          DO J = 1, N
             XX = X(I) + R1 * Y(I)
             YY = R2 * X(I) + Y(I)
             X(I) = XX * B(1,1)
             Y(I) = YY * B(2,2)
             I = I + LDA
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(2,1))
          R2 = B(2,2) / B(2,1)
          !DIR$ VECTOR ALWAYS
          DO J = 1, N
             XX = X(I) + R1 * Y(I)
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
             XX = R1 * X(I) + Y(I)
             YY = R2 * X(I) + Y(I)
             X(I) = XX * B(1,2)
             Y(I) = YY * B(2,2)
             I = I + LDA
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(2,1))
          R2 = B(2,2) / B(2,1)
          !DIR$ VECTOR ALWAYS
          DO J = 1, N
             XX = R1 * X(I) + Y(I)
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
    COMPLEX(KIND=DWP), INTENT(IN) :: B(2,2)
    INTEGER, INTENT(IN) :: M
    COMPLEX(KIND=DWP), INTENT(INOUT) :: X(M), Y(M)

    COMPLEX(KIND=DWP) :: XX, YY
    INTEGER :: I

    !DIR$ VECTOR ALWAYS ASSERT
    DO I = 1, M
       XX = X(I) + Y(I) * B(2,1)
       YY = X(I) * B(1,2) + Y(I)
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
       R1 = B(2,1) / B(1,1)
       IF (ABS(B(2,2)) .GE. ABS(B(1,2))) THEN
          R2 = B(1,2) / B(2,2)
          !DIR$ VECTOR ALWAYS ASSERT
          DO I = 1, M
             XX = X(I) + Y(I) * R1
             YY = X(I) * R2 + Y(I)
             X(I) = XX * B(1,1)
             Y(I) = YY * B(2,2)
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(1,2))
          R2 = B(2,2) / B(1,2)
          !DIR$ VECTOR ALWAYS ASSERT
          DO I = 1, M
             XX = X(I) + Y(I) * R1
             YY = X(I) + Y(I) * R2
             X(I) = XX * B(1,1)
             Y(I) = YY * B(1,2)
          END DO
       END IF
    ELSE ! ABS(B(1,1)) .LT. ABS(B(2,1))
       R1 = B(1,1) / B(2,1)
       IF (ABS(B(2,2)) .GE. ABS(B(1,2))) THEN
          R2 = B(1,2) / B(2,2)
          !DIR$ VECTOR ALWAYS ASSERT
          DO I = 1, M
             XX = X(I) * R1 + Y(I)
             YY = X(I) * R2 + Y(I)
             X(I) = XX * B(2,1)
             Y(I) = YY * B(2,2)
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(1,2))
          R2 = B(2,2) / B(1,2)
          !DIR$ VECTOR ALWAYS ASSERT
          DO I = 1, M
             XX = X(I) * R1 + Y(I)
             YY = X(I) + Y(I) * R2
             X(I) = XX * B(2,1)
             Y(I) = YY * B(1,2)
          END DO
       END IF
    END IF
  END SUBROUTINE AB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE REAL(KIND=DWP) FUNCTION ABODNZ(B, M, X, Y, P, Q)
    ! skipping the diagonal elements X(P) and Y(Q),
    ! D = SUM(|oldX(I)|-|newX(I)| + |oldY(I)|-|newY(I)|)
    !   = ||oldX oldY||_1 - ||newX newY||_1
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: M, P, Q
    COMPLEX(KIND=DWP), INTENT(IN) :: B(2,2), X(M), Y(M)

    COMPLEX(KIND=DWP) :: R1, R2, XX, YY
    INTEGER :: I

    IF ((M .LT. 2) .OR. (P .LE. 0) .OR. (P .GT. M) .OR. (Q .LE. P) .OR. (Q .GT. M)) THEN
       ABODNZ = D_MZERO
       RETURN
    END IF
    ABODNZ = ABS(Y(P)) + ABS(X(Q))

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
             ABODNZ = ABODNZ + (ABS(X(I)) - ABS(XX)) + (ABS(Y(I)) - ABS(YY))
          END DO
          !DIR$ VECTOR ALWAYS
          DO I = P+1, Q-1
             !DIR$ FMA
             XX = (X(I) + Y(I) * R1) * B(1,1)
             !DIR$ FMA
             YY = (X(I) * R2 + Y(I)) * B(2,2)
             ABODNZ = ABODNZ + (ABS(X(I)) - ABS(XX)) + (ABS(Y(I)) - ABS(YY))
          END DO
          !DIR$ VECTOR ALWAYS
          DO I = Q+1, M
             !DIR$ FMA
             XX = (X(I) + Y(I) * R1) * B(1,1)
             !DIR$ FMA
             YY = (X(I) * R2 + Y(I)) * B(2,2)
             ABODNZ = ABODNZ + (ABS(X(I)) - ABS(XX)) + (ABS(Y(I)) - ABS(YY))
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(1,2))
          R2 = B(2,2) / B(1,2)
          !DIR$ VECTOR ALWAYS
          DO I = 1, P-1
             !DIR$ FMA
             XX = (X(I) + Y(I) * R1) * B(1,1)
             !DIR$ FMA
             YY = (X(I) + Y(I) * R2) * B(1,2)
             ABODNZ = ABODNZ + (ABS(X(I)) - ABS(XX)) + (ABS(Y(I)) - ABS(YY))
          END DO
          !DIR$ VECTOR ALWAYS
          DO I = P+1, Q-1
             !DIR$ FMA
             XX = (X(I) + Y(I) * R1) * B(1,1)
             !DIR$ FMA
             YY = (X(I) + Y(I) * R2) * B(1,2)
             ABODNZ = ABODNZ + (ABS(X(I)) - ABS(XX)) + (ABS(Y(I)) - ABS(YY))
          END DO
          !DIR$ VECTOR ALWAYS
          DO I = Q+1, M
             !DIR$ FMA
             XX = (X(I) + Y(I) * R1) * B(1,1)
             !DIR$ FMA
             YY = (X(I) + Y(I) * R2) * B(1,2)
             ABODNZ = ABODNZ + (ABS(X(I)) - ABS(XX)) + (ABS(Y(I)) - ABS(YY))
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
             ABODNZ = ABODNZ + (ABS(X(I)) - ABS(XX)) + (ABS(Y(I)) - ABS(YY))
          END DO
          !DIR$ VECTOR ALWAYS
          DO I = P+1, Q-1
             !DIR$ FMA
             XX = (X(I) * R1 + Y(I)) * B(2,1)
             !DIR$ FMA
             YY = (X(I) * R2 + Y(I)) * B(2,2)
             ABODNZ = ABODNZ + (ABS(X(I)) - ABS(XX)) + (ABS(Y(I)) - ABS(YY))
          END DO
          !DIR$ VECTOR ALWAYS
          DO I = Q+1, M
             !DIR$ FMA
             XX = (X(I) * R1 + Y(I)) * B(2,1)
             !DIR$ FMA
             YY = (X(I) * R2 + Y(I)) * B(2,2)
             ABODNZ = ABODNZ + (ABS(X(I)) - ABS(XX)) + (ABS(Y(I)) - ABS(YY))
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(1,2))
          R2 = B(2,2) / B(1,2)
          !DIR$ VECTOR ALWAYS
          DO I = 1, P-1
             !DIR$ FMA
             XX = (X(I) * R1 + Y(I)) * B(2,1)
             !DIR$ FMA
             YY = (X(I) + Y(I) * R2) * B(1,2)
             ABODNZ = ABODNZ + (ABS(X(I)) - ABS(XX)) + (ABS(Y(I)) - ABS(YY))
          END DO
          !DIR$ VECTOR ALWAYS
          DO I = P+1, Q-1
             !DIR$ FMA
             XX = (X(I) * R1 + Y(I)) * B(2,1)
             !DIR$ FMA
             YY = (X(I) + Y(I) * R2) * B(1,2)
             ABODNZ = ABODNZ + (ABS(X(I)) - ABS(XX)) + (ABS(Y(I)) - ABS(YY))
          END DO
          !DIR$ VECTOR ALWAYS
          DO I = Q+1, M
             !DIR$ FMA
             XX = (X(I) * R1 + Y(I)) * B(2,1)
             !DIR$ FMA
             YY = (X(I) + Y(I) * R2) * B(1,2)
             ABODNZ = ABODNZ + (ABS(X(I)) - ABS(XX)) + (ABS(Y(I)) - ABS(YY))
          END DO
       END IF
    END IF
  END FUNCTION ABODNZ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE ZHSVD2D(A, U, INFO)
    ! A diagonal
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    COMPLEX(KIND=DWP) :: V
    REAL(KIND=DWP) :: W

    INFO = 0

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
       U(1,1) = V * U(1,1)
       U(1,2) = V * U(1,2)
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
       U(2,1) = V * U(2,1)
       U(2,2) = V * U(2,2)
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

    COMPLEX(KIND=DWP) :: V(2,2), W(2,2)

    INFO = 0

    ! TODO: transform
    IF (H) THEN
       CONTINUE
    ELSE
       CONTINUE
    END IF

    CALL CA(V, 2, U(1,1), U(2,1), 2)
    CALL CA(V, 2, A(1,1), A(2,1), 2)

    CALL AC(W, 2, A(1,1), A(1,2))
    CALL AC(W, 2, Z(1,1), Z(1,2))

    CALL ZHSVD2D(A, U, INFO)
  END SUBROUTINE ZHSVD2U

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE ZHSVD2L(A, U, Z, INFO)
    ! A lower triangular, not diagonal, hyperbolic J
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    COMPLEX(KIND=DWP) :: V(2,2), W(2,2)

    INFO = 0

    ! TODO: transform

    CALL CA(V, 2, U(1,1), U(2,1), 2)
    CALL CA(V, 2, A(1,1), A(2,1), 2)

    CALL AC(W, 2, A(1,1), A(1,2))
    CALL AC(W, 2, Z(1,1), Z(1,2))

    CALL ZHSVD2D(A, U, INFO)
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
    S = A(2,1) / A(1,1)
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
       U(1,1) = R * U(1,1)
       U(1,2) = R * U(1,2)
       A(1,2) = R * A(1,2)
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
       U(2,1) = R * U(2,1)
       U(2,2) = R * U(2,2)
    END IF

    IF (A(1,2) .EQ. Z_ZERO) THEN
       ! A diagonal
       IF (H .AND. (INFO .EQ. 1)) THEN
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
          ! early exit
          INFO = 0
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

    ! record in INFO if old diag(A) and new diag(A) differ
    IF ((A(1,1) .NE. A(2,1)) .OR. (A(2,2) .NE. A(1,2))) THEN
       INFO = 1
    ELSE ! old diag(A) .EQ. new diag(A)
       INFO = 0
    END IF

    ! explicitly make A diagonal
    A(2,1) = Z_ZERO
    A(1,2) = Z_ZERO

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

    ! store diag(A) to W
    I = -K
    W(1,1) = SCALE(REAL(A(1,1)), I)
    W(2,1) = SCALE(AIMAG(A(1,1)), I)
    W(1,2) = SCALE(REAL(A(2,2)), I)
    W(2,2) = SCALE(AIMAG(A(2,2)), I)

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
    A(2,1) = CMPLX(W(1,1), W(2,1), DWP)
    A(1,2) = CMPLX(W(1,2), W(2,2), DWP)

    CALL ZHSVD2S((J(1) + J(2)), A, U, Z, INFO)
  END SUBROUTINE ZHSVD2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE ZTRANSF
