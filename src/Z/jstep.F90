MODULE JSTEP
  USE DTYPES
  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE FUNCTION PQI(P, Q)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: P, Q
    INTEGER :: PQI

    IF ((Q .LT. 2) .OR. (P .LT. 1) .OR. (P .GE. Q)) THEN
       PQI = 0
    ELSE
       PQI = ((Q - 2) * (Q - 1)) / 2 + P
    END IF
  END FUNCTION PQI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE FUNCTION AMAG1(APP, AQP, APQ, AQQ, JP, JQ)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(IN) :: APP, AQP, APQ, AQQ
    INTEGER, INTENT(IN) :: JP, JQ
    REAL(KIND=DWP) :: AMAG1

    AMAG1 = ABS(AQP) + ABS(APQ)
  END FUNCTION AMAG1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE BUILD_STEP(N, A, LDA, J, AM, DZ, S, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, LDA, J(N)
    COMPLEX(KIND=DWP), INTENT(IN) :: A(LDA,N)
    PROCEDURE(AMAG) :: AM
    TYPE(DZBW), INTENT(OUT), TARGET :: DZ((N*(N-1))/2)
    INTEGER, INTENT(OUT) :: S(N/2), INFO

    INTEGER :: NN, N_2, P, Q, I

    IF (N .LT. 0) THEN
       INFO = -1
    ELSE IF (LDA .LT. N) THEN
       INFO = -3
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN
    IF (N .EQ. 0) RETURN

    I = 1
    DO Q = 2, N
       DO P = 1, Q-1
          DZ(I)%P = P
          DZ(I)%Q = Q
          I = I + 1
       END DO
    END DO

    NN = (N * (N - 1)) / 2
    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(P,Q,I) SHARED(N,A,J,DZ)
    DO I = 1, NN
       P = DZ(I)%P
       Q = DZ(I)%Q
       CALL DZBW_GEN(A(P,P), A(Q,P), A(P,Q), A(Q,Q), J(P), J(Q), P, Q, AM, DZ(I))
    END DO
    !$OMP END PARALLEL DO

    CALL DZBW_SRT(NN, DZ, INFO)
    IF (INFO .NE. 0) RETURN

    N_2 = N / 2
    CALL DZBW_NCP(NN, DZ, N_2, S, INFO)
  END SUBROUTINE BUILD_STEP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE JSTEP
