!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE REAL(KIND=DWP) FUNCTION ZMAG1(N, P, Q, A, LDA, J)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, P, Q, LDA, J(N)
    COMPLEX(KIND=DWP), INTENT(IN) :: A(LDA,N)

    IF ((A(Q,P) .NE. Z_ZERO) .OR. (A(P,Q) .NE. Z_ZERO) .OR. (AIMAG(A(P,P)) .NE. D_ZERO) .OR. (AIMAG(A(Q,Q)) .NE. D_ZERO) .OR. &
         (SIGN(D_ONE, REAL(A(P,P))) .EQ. D_MONE) .OR. (SIGN(D_ONE, REAL(A(Q,Q))) .EQ. D_MONE) .OR. &
         ((J(P) .EQ. J(Q)) .AND. (REAL(A(P,P)) .LT. REAL(A(Q,Q))))) THEN
       ZMAG1 = D_ZERO
       !DIR$ FMA
       ZMAG1 = ZMAG1 + REAL(A(Q,P)) * REAL(A(Q,P))
       !DIR$ FMA
       ZMAG1 = ZMAG1 + AIMAG(A(Q,P)) * AIMAG(A(Q,P))
       !DIR$ FMA
       ZMAG1 = ZMAG1 + REAL(A(P,Q)) * REAL(A(P,Q))
       !DIR$ FMA
       ZMAG1 = ZMAG1 + AIMAG(A(P,Q)) * AIMAG(A(P,Q))
       ! not to be used in this form, but only as a mark of a hyperbolic case
       ! IF (J(P) .NE. J(Q)) ZMAG1 = -ZMAG1
    ELSE ! no transform
       ZMAG1 = QUIET_NAN((P - 1) * N + (Q - 1))
    END IF
  END FUNCTION ZMAG1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE REAL(KIND=DWP) FUNCTION DMAG1(N, P, Q, A, LDA, J)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, P, Q, LDA, J(N)
    REAL(KIND=DWP), INTENT(IN) :: A(LDA,N)

    REAL(KIND=DWP) :: A2(2,2), U2(2,2), Z2(2,2)
    INTEGER :: INFO

    IF ((A(Q,P) .NE. D_ZERO) .OR. (A(P,Q) .NE. D_ZERO) .OR. (SIGN(D_ONE, A(P,P)) .EQ. D_MONE) .OR. &
         (SIGN(D_ONE, A(Q,Q)) .EQ. D_MONE) .OR. ((J(P) .EQ. J(Q)) .AND. (A(P,P) .LT. A(Q,Q)))) THEN
       IF (J(P) .EQ. J(Q)) THEN
          DMAG1 = D_ZERO
          !DIR$ FMA
          DMAG1 = DMAG1 + A(Q,P) * A(Q,P)
          !DIR$ FMA
          DMAG1 = DMAG1 + A(P,Q) * A(P,Q)
       ELSE ! J(P) .NE. J(Q)
          A2(1,1) = A(P,P)
          A2(2,1) = A(Q,P)
          A2(1,2) = A(P,Q)
          A2(2,2) = A(Q,Q)
          CALL DHSVD2(.TRUE., A2, U2, Z2, INFO)
          IF (INFO .LE. 0) THEN
             DMAG1 = QUIET_NAN((P - 1) * N + (Q - 1))
          ELSE ! a non-trivial transform
             DMAG1 = ABODND(Z2, N, A(1,P), A(1,Q), P, Q)
          END IF
       END IF
    ELSE ! no transform
       DMAG1 = QUIET_NAN((P - 1) * N + (Q - 1))
    END IF
  END FUNCTION DMAG1

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
