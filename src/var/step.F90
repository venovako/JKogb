  REAL(KIND=DWP), PARAMETER :: DZEPS = D_EPS !* SQRT(D_TWO)
  REAL(KIND=DWP), PARAMETER :: D_RSQRT2 = SCALE(SQRT(D_TWO), -1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE REAL(KIND=DWP) FUNCTION ZMAG2(N, P, Q, A, LDA, J)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, P, Q, LDA, J(N)
    COMPLEX(KIND=DWP), INTENT(IN) :: A(LDA,N)

    REAL(KIND=DWP) :: AAPP, AAQP, AAPQ, AAQQ, MAXPQ, MINPQ

    AAPP = ABS(A(P,P))
    AAQP = ABS(A(Q,P))
    AAPQ = ABS(A(P,Q))
    AAQQ = ABS(A(Q,Q))

    IF (AAPP .GE. AAQQ) THEN
       MAXPQ = AAPP
       MINPQ = AAQQ
    ELSE
       MAXPQ = AAQQ
       MINPQ = AAPP
    END IF

    IF (MAX(AAQP, AAPQ) .GT. ((MAXPQ * DZEPS) * MINPQ)) THEN
       ZMAG2 = AAQP + AAPQ
    ELSE ! no transform
       ZMAG2 = QUIET_NAN((P - 1) * N + (Q - 1))
    END IF
  END FUNCTION ZMAG2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE REAL(KIND=DWP) FUNCTION DMAG2(N, P, Q, A, LDA, J)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, P, Q, LDA, J(N)
    REAL(KIND=DWP), INTENT(IN) :: A(LDA,N)

    REAL(KIND=DWP) :: AAPP, AAQP, AAPQ, AAQQ, MAXPQ, MINPQ

    AAPP = ABS(A(P,P))
    AAQP = ABS(A(Q,P))
    AAPQ = ABS(A(P,Q))
    AAQQ = ABS(A(Q,Q))

    IF (AAPP .GE. AAQQ) THEN
       MAXPQ = AAPP
       MINPQ = AAQQ
    ELSE
       MAXPQ = AAQQ
       MINPQ = AAPP
    END IF

    IF (MAX(AAQP, AAPQ) .GT. ((MAXPQ * D_EPS) * MINPQ)) THEN
       DMAG2 = AAQP + AAPQ
    ELSE ! no transform
       DMAG2 = QUIET_NAN((P - 1) * N + (Q - 1))
    END IF
  END FUNCTION DMAG2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
