MODULE dtransf
  USE ktransf
#ifndef FMAD
#define FMAD(a,b,c) ((a)*(b)+(c))
#endif
  IMPLICIT NONE

CONTAINS

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

  PURE SUBROUTINE UT2(U)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(INOUT) :: U(2,2)

    REAL(KIND=DWP) :: U21

    U21 = U(2,1)
    U(2,1) = U(1,2)
    U(1,2) = U21
  END SUBROUTINE UT2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE JZTJ2(Z, J)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(INOUT) :: Z(2,2)
    INTEGER, INTENT(IN) :: J(2)

    CALL UT2(Z)
    IF (J(1) .EQ. -1) THEN
       Z(2,1) = -Z(2,1)
       Z(1,2) = -Z(1,2)
    END IF
    IF (J(2) .EQ. -1) THEN
       Z(2,1) = -Z(2,1)
       Z(1,2) = -Z(1,2)
    END IF
  END SUBROUTINE JZTJ2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE C1A(B, A)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: B(2,2)
    REAL(KIND=DWP), INTENT(INOUT) :: A(2)

    REAL(KIND=DWP) :: C(2)

    C(1) = FMAD(B(1,2), A(2), A(1)) / B(1,1)
    C(2) = FMAD(B(2,1), A(1), A(2)) / B(2,2)
    A = C
  END SUBROUTINE C1A

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE C2A(B, A)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: B(2,2)
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2)

    REAL(KIND=DWP) :: C(2,2)
    INTEGER :: J

    DO J = 1, 2
       C(1,J) = FMAD(B(1,2), A(2,J), A(1,J)) / B(1,1)
       C(2,J) = FMAD(B(2,1), A(1,J), A(2,J)) / B(2,2)
    END DO
    A = C
  END SUBROUTINE C2A

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE A2C(A, B)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2)
    REAL(KIND=DWP), INTENT(IN) :: B(2,2)

    REAL(KIND=DWP) :: C(2,2)
    INTEGER :: I

    DO I = 1, 2
       C(I,1) = FMAD(A(I,2), B(2,1), A(I,1)) / B(1,1)
       C(I,2) = FMAD(A(I,1), B(1,2), A(I,2)) / B(2,2)
    END DO
    A = C
  END SUBROUTINE A2C

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
          DO J = 1, N
             XX = FMAD(R1, Y(I), X(I))
             YY = FMAD(R2, X(I), Y(I))
             X(I) = XX * B(1,1)
             Y(I) = YY * B(2,2)
             I = I + LDA
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(2,1))
          R2 = B(2,2) / B(2,1)
          DO J = 1, N
             XX = FMAD(R1, Y(I), X(I))
             YY = FMAD(R2, Y(I), X(I))
             X(I) = XX * B(1,1)
             Y(I) = YY * B(2,1)
             I = I + LDA
          END DO
       END IF
    ELSE ! ABS(B(1,1)) .LT. ABS(B(1,2))
       R1 = B(1,1) / B(1,2)
       IF (ABS(B(2,2)) .GE. ABS(B(2,1))) THEN
          R2 = B(2,1) / B(2,2)
          DO J = 1, N
             XX = FMAD(R1, X(I), Y(I))
             YY = FMAD(R2, X(I), Y(I))
             X(I) = XX * B(1,2)
             Y(I) = YY * B(2,2)
             I = I + LDA
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(2,1))
          R2 = B(2,2) / B(2,1)
          DO J = 1, N
             XX = FMAD(R1, X(I), Y(I))
             YY = FMAD(R2, Y(I), X(I))
             X(I) = XX * B(1,2)
             Y(I) = YY * B(2,1)
             I = I + LDA
          END DO
       END IF
    END IF
  END SUBROUTINE BA

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
          DO I = 1, M
             XX = FMAD(Y(I), R1, X(I))
             YY = FMAD(X(I), R2, Y(I))
             X(I) = XX * B(1,1)
             Y(I) = YY * B(2,2)
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(1,2))
          R2 = B(2,2) / B(1,2)
          DO I = 1, M
             XX = FMAD(Y(I), R1, X(I))
             YY = FMAD(Y(I), R2, X(I))
             X(I) = XX * B(1,1)
             Y(I) = YY * B(1,2)
          END DO
       END IF
    ELSE ! ABS(B(1,1)) .LT. ABS(B(2,1))
       R1 = B(1,1) / B(2,1)
       IF (ABS(B(2,2)) .GE. ABS(B(1,2))) THEN
          R2 = B(1,2) / B(2,2)
          DO I = 1, M
             XX = FMAD(X(I), R1, Y(I))
             YY = FMAD(X(I), R2, Y(I))
             X(I) = XX * B(2,1)
             Y(I) = YY * B(2,2)
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(1,2))
          R2 = B(2,2) / B(1,2)
          DO I = 1, M
             XX = FMAD(X(I), R1, Y(I))
             YY = FMAD(Y(I), R2, X(I))
             X(I) = XX * B(2,1)
             Y(I) = YY * B(1,2)
          END DO
       END IF
    END IF
  END SUBROUTINE AB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE REAL(KIND=DWP) FUNCTION ABODNDF2(B, M, X, Y, P, Q)
    ! skipping the diagonal elements X(P) and Y(Q),
    ! D = SUM(oldX(I)^2-newX(I)^2 + oldY(I)^2-newY(I)^2)
    !   = ||oldX oldY||_F^2 - ||newX newY||_F^2
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: M, P, Q
    REAL(KIND=DWP), INTENT(IN) :: B(2,2), X(M), Y(M)

    REAL(KIND=DWP) :: R1, R2, XX, YY
    INTEGER :: I

    ABODNDF2 = DASUM2(X(Q), Y(P))

    IF (ABS(B(1,1)) .GE. ABS(B(2,1))) THEN
       R1 = B(2,1) / B(1,1)
       IF (ABS(B(2,2)) .GE. ABS(B(1,2))) THEN
          R2 = B(1,2) / B(2,2)
          DO I = 1, P-1
             XX = FMAD(Y(I), R1, X(I)) * B(1,1)
             YY = FMAD(X(I), R2, Y(I)) * B(2,2)
             ABODNDF2 = FMAD((X(I) - XX), (X(I) + XX), ABODNDF2)
             ABODNDF2 = FMAD((Y(I) - YY), (Y(I) + YY), ABODNDF2)
          END DO
          DO I = P+1, Q-1
             XX = FMAD(Y(I), R1, X(I)) * B(1,1)
             YY = FMAD(X(I), R2, Y(I)) * B(2,2)
             ABODNDF2 = FMAD((X(I) - XX), (X(I) + XX), ABODNDF2)
             ABODNDF2 = FMAD((Y(I) - YY), (Y(I) + YY), ABODNDF2)
          END DO
          DO I = Q+1, M
             XX = FMAD(Y(I), R1, X(I)) * B(1,1)
             YY = FMAD(X(I), R2, Y(I)) * B(2,2)
             ABODNDF2 = FMAD((X(I) - XX), (X(I) + XX), ABODNDF2)
             ABODNDF2 = FMAD((Y(I) - YY), (Y(I) + YY), ABODNDF2)
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(1,2))
          R2 = B(2,2) / B(1,2)
          DO I = 1, P-1
             XX = FMAD(Y(I), R1, X(I)) * B(1,1)
             YY = FMAD(Y(I), R2, X(I)) * B(1,2)
             ABODNDF2 = FMAD((X(I) - XX), (X(I) + XX), ABODNDF2)
             ABODNDF2 = FMAD((Y(I) - YY), (Y(I) + YY), ABODNDF2)
          END DO
          DO I = P+1, Q-1
             XX = FMAD(Y(I), R1, X(I)) * B(1,1)
             YY = FMAD(Y(I), R2, X(I)) * B(1,2)
             ABODNDF2 = FMAD((X(I) - XX), (X(I) + XX), ABODNDF2)
             ABODNDF2 = FMAD((Y(I) - YY), (Y(I) + YY), ABODNDF2)
          END DO
          DO I = Q+1, M
             XX = FMAD(Y(I), R1, X(I)) * B(1,1)
             YY = FMAD(Y(I), R2, X(I)) * B(1,2)
             ABODNDF2 = FMAD((X(I) - XX), (X(I) + XX), ABODNDF2)
             ABODNDF2 = FMAD((Y(I) - YY), (Y(I) + YY), ABODNDF2)
          END DO
       END IF
    ELSE ! ABS(B(1,1)) .LT. ABS(B(2,1))
       R1 = B(1,1) / B(2,1)
       IF (ABS(B(2,2)) .GE. ABS(B(1,2))) THEN
          R2 = B(1,2) / B(2,2)
          DO I = 1, P-1
             XX = FMAD(X(I), R1, Y(I)) * B(2,1)
             YY = FMAD(X(I), R2, Y(I)) * B(2,2)
             ABODNDF2 = FMAD((X(I) - XX), (X(I) + XX), ABODNDF2)
             ABODNDF2 = FMAD((Y(I) - YY), (Y(I) + YY), ABODNDF2)
          END DO
          DO I = P+1, Q-1
             XX = FMAD(X(I), R1, Y(I)) * B(2,1)
             YY = FMAD(X(I), R2, Y(I)) * B(2,2)
             ABODNDF2 = FMAD((X(I) - XX), (X(I) + XX), ABODNDF2)
             ABODNDF2 = FMAD((Y(I) - YY), (Y(I) + YY), ABODNDF2)
          END DO
          DO I = Q+1, M
             XX = FMAD(X(I), R1, Y(I)) * B(2,1)
             YY = FMAD(X(I), R2, Y(I)) * B(2,2)
             ABODNDF2 = FMAD((X(I) - XX), (X(I) + XX), ABODNDF2)
             ABODNDF2 = FMAD((Y(I) - YY), (Y(I) + YY), ABODNDF2)
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(1,2))
          R2 = B(2,2) / B(1,2)
          DO I = 1, P-1
             XX = FMAD(X(I), R1, Y(I)) * B(2,1)
             YY = FMAD(Y(I), R2, X(I)) * B(1,2)
             ABODNDF2 = FMAD((X(I) - XX), (X(I) + XX), ABODNDF2)
             ABODNDF2 = FMAD((Y(I) - YY), (Y(I) + YY), ABODNDF2)
          END DO
          DO I = P+1, Q-1
             XX = FMAD(X(I), R1, Y(I)) * B(2,1)
             YY = FMAD(Y(I), R2, X(I)) * B(1,2)
             ABODNDF2 = FMAD((X(I) - XX), (X(I) + XX), ABODNDF2)
             ABODNDF2 = FMAD((Y(I) - YY), (Y(I) + YY), ABODNDF2)
          END DO
          DO I = Q+1, M
             XX = FMAD(X(I), R1, Y(I)) * B(2,1)
             YY = FMAD(Y(I), R2, X(I)) * B(1,2)
             ABODNDF2 = FMAD((X(I) - XX), (X(I) + XX), ABODNDF2)
             ABODNDF2 = FMAD((Y(I) - YY), (Y(I) + YY), ABODNDF2)
          END DO
       END IF
    END IF
  END FUNCTION ABODNDF2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE DHSVD2T(A, J, U, Z, INFO)
    ! A general
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(IN) :: J(2)
    INTEGER, INTENT(OUT) :: INFO

    REAL(KIND=DWP) :: Q(2,2), C, S, R
    LOGICAL :: P, D

    D = .TRUE.
    IF (A(2,1) .EQ. D_ZERO) THEN
       C = ABS(A(1,1))
    ELSE IF (A(1,1) .EQ. D_ZERO) THEN
       C = ABS(A(2,1))
       D = .FALSE.
    ELSE ! full 1st column
       C = HYPOT(A(1,1), A(2,1))
       D = .FALSE.
    END IF
    IF (A(1,2) .EQ. D_ZERO) THEN
       S = ABS(A(2,2))
    ELSE IF (A(2,2) .EQ. D_ZERO) THEN
       S = ABS(A(1,2))
       D = .FALSE.
    ELSE ! full 2nd column
       S = HYPOT(A(1,2), A(2,2))
       D = .FALSE.
    END IF

    ! column pivoting
    IF (C .LT. S) THEN
       R = S
       ! swap the columns of A
       C = A(1,1)
       S = A(2,1)
       A(1,1) = A(1,2)
       A(2,1) = A(2,2)
       A(1,2) = C
       A(2,2) = S
       ! swap the columns of Z
       IF (J(1) .EQ. J(2)) THEN
          C = Z(1,1)
          S = Z(2,1)
          Z(1,1) = Z(1,2)
          Z(2,1) = Z(2,2)
          Z(1,2) = C
          Z(2,2) = S
       END IF
       ! record the swap
       P = .TRUE.
    ELSE ! no pivoting
       R = C
       P = .FALSE.
    END IF

    ! make A(1,1) non-negative
    IF (SIGN(D_ONE, A(1,1)) .EQ. D_MONE) THEN
       U(1,1) = -U(1,1)
       U(1,2) = -U(1,2)
       A(1,2) = -A(1,2)
       A(1,1) = -A(1,1)
    END IF

    ! make A(2,1) non-negative
    IF (SIGN(D_ONE, A(2,1)) .EQ. D_MONE) THEN
       U(2,1) = -U(2,1)
       U(2,2) = -U(2,2)
       A(2,2) = -A(2,2)
       A(2,1) = -A(2,1)
    END IF

    ! row sorting
    IF (A(1,1) .LT. A(2,1)) THEN
       ! swap the rows of U
       C = U(1,1)
       S = U(1,2)
       U(1,1) = U(2,1)
       U(1,2) = U(2,2)
       U(2,1) = C
       U(2,2) = S
       ! swap the rows of A
       C = A(1,1)
       S = A(1,2)
       A(1,1) = A(2,1)
       A(1,2) = A(2,2)
       A(2,1) = C
       A(2,2) = S
    END IF

    ! QR factorization of A
    ! |A(1,1)| >= |A(2,1)| cannot be 0; otherwise,
    ! ||A_1||=0 >= ||A_2|| >= 0, so A=0
    ! specifically, A is diagonal

    ! S is tangent here, |S| <= 1
    IF (.NOT. D) THEN
       S = A(2,1) / A(1,1)
       C = SQRT(FMAD(S, S, D_ONE)) ! D_ONE /
       Q(1,1) =  C
       Q(2,1) = -S
       Q(1,2) =  S
       Q(2,2) =  C

       ! premultiply U by Q^T
       CALL C2A(Q, U)
       CALL C1A(Q, A(1,2))
    END IF
    A(1,1) = R
    A(2,1) = D_ZERO

    IF ((J(1) .NE. J(2)) .AND. P) THEN
       ! swap the columns of A
       C = A(1,1)
       S = A(2,1)
       A(1,1) = A(1,2)
       A(2,1) = A(2,2)
       A(1,2) = C
       A(2,2) = S
       ! A upper antitriangular
       !     | x R | <- R is the largest
       ! A = | y 0 |    element by magnitude
       ! swap the rows of U
       C = U(1,1)
       S = U(1,2)
       U(1,1) = U(2,1)
       U(1,2) = U(2,2)
       U(2,1) = C
       U(2,2) = S
       ! swap the rows of A
       C = A(1,1)
       S = A(1,2)
       A(1,1) = A(2,1)
       A(1,2) = A(2,2)
       A(2,1) = C
       A(2,2) = S
       ! A lower triangular
       ! A = | y 0 |    R is the largest
       !     | x R | <- element by magnitude

       ! make A(2,1) non-negative
       IF (SIGN(D_ONE, A(2,1)) .EQ. D_MONE) THEN
          Z(1,1) = -Z(1,1)
          Z(2,1) = -Z(2,1)
          A(1,1) = -A(1,1)
          A(2,1) = -A(2,1)
       END IF

       ! make A(1,1) non-negative
       IF (SIGN(D_ONE, A(1,1)) .EQ. D_MONE) THEN
          U(1,1) = -U(1,1)
          U(1,2) = -U(1,2)
          A(1,1) = -A(1,1)
       END IF
    ELSE ! trigonometric or no column pivoting
       ! A upper triangular
       !     | R x | <- R is the largest
       ! A = | 0 y |    element by magnitude

       ! make A(1,2) non-negative
       IF (SIGN(D_ONE, A(1,2)) .EQ. D_MONE) THEN
          Z(1,2) = -Z(1,2)
          Z(2,2) = -Z(2,2)
          A(2,2) = -A(2,2)
          A(1,2) = -A(1,2)
       END IF

       ! make A(2,2) non-negative
       IF (SIGN(D_ONE, A(2,2)) .EQ. D_MONE) THEN
          U(2,1) = -U(2,1)
          U(2,2) = -U(2,2)
          A(2,2) = -A(2,2)
       END IF
    END IF

    IF (D) THEN
       INFO = 0
    ELSE ! .NOT. D
       INFO = 1
    END IF
  END SUBROUTINE DHSVD2T

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE DHSVD2S(A, J, U, Z, INFO)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(IN) :: J(2)
    INTEGER, INTENT(INOUT) :: INFO

    INTEGER :: T

    T = J(1) + J(2)

    IF (((T .EQ. 2) .AND. (A(1,1) .LT. A(2,2))) .OR. ((T .EQ. -2) .AND. (A(2,2) .LT. A(1,1)))) THEN
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
    A(2,1) = D_ZERO
    A(1,2) = D_ZERO

    ! check if U is identity and record in INFO if it is not
    IF ((U(1,1) .NE. D_ONE) .OR. (U(2,1) .NE. D_ZERO) .OR. (U(1,2) .NE. D_ZERO) .OR. (U(2,2) .NE. D_ONE)) INFO = INFO + 2
    ! check if Z is identity and record in INFO if it is not
    IF ((Z(1,1) .NE. D_ONE) .OR. (Z(2,1) .NE. D_ZERO) .OR. (Z(1,2) .NE. D_ZERO) .OR. (Z(2,2) .NE. D_ONE)) INFO = INFO + 4
    ! check for overflow
    ! IF (A(1,1) .GT. HUGE(D_ZERO)) INFO = INFO + 8
    ! IF (A(2,2) .GT. HUGE(D_ZERO)) INFO = INFO + 16
  END SUBROUTINE DHSVD2S

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE DSCALEA(A, S)
    IMPLICIT NONE

    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2)
    INTEGER, INTENT(OUT) :: S

    REAL(KIND=DWP) :: AA
    INTEGER :: I, J

    S = HUGE(S)
    DO J = 1, 2
       DO I = 1, 2
          AA = ABS(A(I,J))
          IF (.NOT. (AA .LE. HUGE(D_ZERO))) THEN
             ! -3 <= S <= -6
             S = -(J * 2 + I)
             RETURN
          END IF
          S = MIN(S, (EH - EXPONENT(MAX(AA,MINF)) - 2))
       END DO
    END DO

    IF (S .NE. 0) THEN
       DO J = 1, 2
          DO I = 1, 2
             A(I,J) = SCALE(A(I,J), S)
          END DO
       END DO
    END IF
  END SUBROUTINE DSCALEA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE DHSVD2(A, J, U, Z, INFO)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2)
    INTEGER, INTENT(IN) :: J(2)
    REAL(KIND=DWP), INTENT(OUT) :: U(2,2), Z(2,2)
    INTEGER, INTENT(OUT) :: INFO

    REAL(KIND=DWP) :: B(2,2), C(2,2)
    INTEGER :: S

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
    ELSE ! A OK
       INFO = 0
#endif
    END IF

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

    CALL DHSVD2T(A, J, U, Z, INFO)
    IF (INFO .GT. 0) THEN
       CALL KHSVD2(A, J, B, C, S, INFO)
    ELSE IF (INFO .EQ. 0) THEN
       A(1,1) = SCALE(A(1,1), -S)
       A(2,2) = SCALE(A(2,2), -S)
    END IF
    IF (INFO .LT. 0) RETURN
    IF (INFO .GT. 0) THEN
       CALL C2A(B, U)
       CALL A2C(Z, C)
    END IF
    CALL DHSVD2S(A, J, U, Z, INFO)
  END SUBROUTINE DHSVD2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE dtransf
