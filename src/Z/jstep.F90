MODULE JSTEP
  USE DTYPES
  IMPLICIT NONE

  INTEGER, PARAMETER :: AMAG_CNT = 1

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

  PURE SUBROUTINE INIT_AMAG(N, N_2, ID, AMS, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, N_2, ID
    TYPE(MAXSAM), INTENT(OUT) :: AMS
    INTEGER, INTENT(OUT) :: INFO

    INTEGER :: MAX_N_2

    IF (N .LT. 0) THEN
       INFO = -1
    ELSE ! N >= 0
       MAX_N_2 = N / 2
       IF (N_2 .LT. 0) THEN
          INFO = -2
       ELSE IF (N_2 .GT. MAX_N_2) THEN
          INFO = -2
       ELSE IF (ID .LT. 0) THEN
          INFO = -3
       ELSE IF (ID .GT. AMAG_CNT) THEN
          INFO = -3
       ELSE
          INFO = 0
       END IF
    END IF
    IF (INFO .NE. 0) RETURN

    IF (N_2 .EQ. 0) THEN
       AMS%MAXS = MAX_N_2
    ELSE
       AMS%MAXS = N_2
    END IF

    IF (ID .EQ. 0) THEN
       INFO = 1
    ELSE
       INFO = ID
    END IF

    SELECT CASE (INFO)
    CASE (1)
       AMS%AM => AMAG1
    CASE DEFAULT
       ! should never happen
       AMS%AM => NULL()
       INFO = -3
    END SELECT
  END SUBROUTINE INIT_AMAG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE INIT_TRIU(N, P, Q, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(OUT) :: P((N*(N-1))/2), Q((N*(N-1))/2), INFO

    INTEGER :: IP, IQ, I

    IF (N .LT. 0) THEN
       INFO = -1
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN
    IF (N .LE. 1) RETURN

    INFO = (N * (N - 1)) / 2

    I = 1
    DO IQ = 2, N
       DO IP = 1, IQ-1
          P(I) = IP
          Q(I) = IQ
          I = I + 1
       END DO
    END DO
  END SUBROUTINE INIT_TRIU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE INIT_JSTEP(N, N_2, ID, AMS, P, Q, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(INOUT) :: N_2, ID
    TYPE(MAXSAM), INTENT(OUT) :: AMS
    INTEGER, INTENT(OUT) :: P((N*(N-1))/2), Q((N*(N-1))/2), INFO

    CALL INIT_AMAG(N, N_2, ID, AMS, INFO)
    IF (INFO .LE. 0) RETURN

    N_2 = AMS%MAXS
    ID = INFO

    CALL INIT_TRIU(N, P, Q, INFO)
  END SUBROUTINE INIT_JSTEP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE BUILD_STEP(N, A, LDA, J, NN, P, Q, AM, DZ, N_2, S, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, LDA, J(N), NN, P(NN), Q(NN), N_2
    COMPLEX(KIND=DWP), INTENT(IN) :: A(LDA,N)
    PROCEDURE(AMAG) :: AM
    TYPE(DZBW), INTENT(OUT), TARGET :: DZ(NN)
    INTEGER, INTENT(OUT) :: S(N_2), INFO

    INTEGER :: IP, IQ, I

    IF (N .LT. 0) THEN
       INFO = -1
    ELSE IF (LDA .LT. N) THEN
       INFO = -3
    ELSE IF (NN .LT. 0) THEN
       INFO = -5
    ELSE IF (NN .GT. ((N * (N - 1)) / 2)) THEN
       INFO = -5
    ELSE IF (N_2 .LT. 0) THEN
       INFO = -10
    ELSE IF (N_2 .GT. (N / 2)) THEN
       INFO = -10
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN
    IF (N .EQ. 0) RETURN
    IF (NN .EQ. 0) RETURN
    IF (N_2 .EQ. 0) RETURN

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(IP,IQ,I) SHARED(NN,A,J,P,Q,DZ)
    DO I = 1, NN
       IP = P(I)
       IQ = Q(I)
       CALL DZBW_GEN(A(IP,IP), A(IQ,IP), A(IP,IQ), A(IQ,IQ), J(IP), J(IQ), IP, IQ, AM, DZ(I))
    END DO
    !$OMP END PARALLEL DO

    CALL DZBW_SRT(NN, DZ, INFO)
    IF (INFO .NE. 0) RETURN

    CALL DZBW_NCP(NN, DZ, N_2, S, INFO)
  END SUBROUTINE BUILD_STEP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE JSTEP
