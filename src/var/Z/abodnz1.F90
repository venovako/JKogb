  PURE REAL(KIND=DWP) FUNCTION ABODNZ1(B, M, X, Y, P, Q)
    ! skipping the diagonal elements X(P) and Y(Q),
    ! Z = SUM(|oldX(I)|-|newX(I)| + |oldY(I)|-|newY(I)|)
    !   = ||oldX oldY||_1 - ||newX newY||_1
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: M, P, Q
    COMPLEX(KIND=DWP), INTENT(IN) :: B(2,2), X(M), Y(M)

    COMPLEX(KIND=DWP) :: R1, R2, XX, YY
    INTEGER :: I

    IF ((M .LT. 2) .OR. (P .LE. 0) .OR. (P .GT. M) .OR. (Q .LE. P) .OR. (Q .GT. M)) THEN
       ABODNZ1 = D_MZERO
       RETURN
    END IF
    ABODNZ1 = ABS(Y(P)) + ABS(X(Q))

    IF (ABS(B(1,1)) .GE. ABS(B(2,1))) THEN
       R1 = B(2,1) / B(1,1)
       IF (ABS(B(2,2)) .GE. ABS(B(1,2))) THEN
          R2 = B(1,2) / B(2,2)
          DO I = 1, P-1
             XX = (X(I) + Y(I) * R1) * B(1,1)
             YY = (X(I) * R2 + Y(I)) * B(2,2)
             ABODNZ1 = ABODNZ1 + (ABS(X(I)) - ABS(XX)) + (ABS(Y(I)) - ABS(YY))
          END DO
          DO I = P+1, Q-1
             XX = (X(I) + Y(I) * R1) * B(1,1)
             YY = (X(I) * R2 + Y(I)) * B(2,2)
             ABODNZ1 = ABODNZ1 + (ABS(X(I)) - ABS(XX)) + (ABS(Y(I)) - ABS(YY))
          END DO
          DO I = Q+1, M
             XX = (X(I) + Y(I) * R1) * B(1,1)
             YY = (X(I) * R2 + Y(I)) * B(2,2)
             ABODNZ1 = ABODNZ1 + (ABS(X(I)) - ABS(XX)) + (ABS(Y(I)) - ABS(YY))
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(1,2))
          R2 = B(2,2) / B(1,2)
          DO I = 1, P-1
             XX = (X(I) + Y(I) * R1) * B(1,1)
             YY = (X(I) + Y(I) * R2) * B(1,2)
             ABODNZ1 = ABODNZ1 + (ABS(X(I)) - ABS(XX)) + (ABS(Y(I)) - ABS(YY))
          END DO
          DO I = P+1, Q-1
             XX = (X(I) + Y(I) * R1) * B(1,1)
             YY = (X(I) + Y(I) * R2) * B(1,2)
             ABODNZ1 = ABODNZ1 + (ABS(X(I)) - ABS(XX)) + (ABS(Y(I)) - ABS(YY))
          END DO
          DO I = Q+1, M
             XX = (X(I) + Y(I) * R1) * B(1,1)
             YY = (X(I) + Y(I) * R2) * B(1,2)
             ABODNZ1 = ABODNZ1 + (ABS(X(I)) - ABS(XX)) + (ABS(Y(I)) - ABS(YY))
          END DO
       END IF
    ELSE ! ABS(B(1,1)) .LT. ABS(B(2,1))
       R1 = B(1,1) / B(2,1)
       IF (ABS(B(2,2)) .GE. ABS(B(1,2))) THEN
          R2 = B(1,2) / B(2,2)
          DO I = 1, P-1
             XX = (X(I) * R1 + Y(I)) * B(2,1)
             YY = (X(I) * R2 + Y(I)) * B(2,2)
             ABODNZ1 = ABODNZ1 + (ABS(X(I)) - ABS(XX)) + (ABS(Y(I)) - ABS(YY))
          END DO
          DO I = P+1, Q-1
             XX = (X(I) * R1 + Y(I)) * B(2,1)
             YY = (X(I) * R2 + Y(I)) * B(2,2)
             ABODNZ1 = ABODNZ1 + (ABS(X(I)) - ABS(XX)) + (ABS(Y(I)) - ABS(YY))
          END DO
          DO I = Q+1, M
             XX = (X(I) * R1 + Y(I)) * B(2,1)
             YY = (X(I) * R2 + Y(I)) * B(2,2)
             ABODNZ1 = ABODNZ1 + (ABS(X(I)) - ABS(XX)) + (ABS(Y(I)) - ABS(YY))
          END DO
       ELSE ! ABS(B(2,2)) .LT. ABS(B(1,2))
          R2 = B(2,2) / B(1,2)
          DO I = 1, P-1
             XX = (X(I) * R1 + Y(I)) * B(2,1)
             YY = (X(I) + Y(I) * R2) * B(1,2)
             ABODNZ1 = ABODNZ1 + (ABS(X(I)) - ABS(XX)) + (ABS(Y(I)) - ABS(YY))
          END DO
          DO I = P+1, Q-1
             XX = (X(I) * R1 + Y(I)) * B(2,1)
             YY = (X(I) + Y(I) * R2) * B(1,2)
             ABODNZ1 = ABODNZ1 + (ABS(X(I)) - ABS(XX)) + (ABS(Y(I)) - ABS(YY))
          END DO
          DO I = Q+1, M
             XX = (X(I) * R1 + Y(I)) * B(2,1)
             YY = (X(I) + Y(I) * R2) * B(1,2)
             ABODNZ1 = ABODNZ1 + (ABS(X(I)) - ABS(XX)) + (ABS(Y(I)) - ABS(YY))
          END DO
       END IF
    END IF
  END FUNCTION ABODNZ1
