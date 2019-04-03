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

  PURE SUBROUTINE AMP_INIT(ID, AM, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ID
    TYPE(AMP), INTENT(OUT) :: AM
    INTEGER, INTENT(OUT) :: INFO

    IF (ID .LT. 0) THEN
       INFO = -1
    ELSE IF (ID .EQ. 0) THEN
       INFO = 1
    ELSE
       INFO = ID
    END IF
    IF (INFO .LE. 0) RETURN

    SELECT CASE (INFO)
    CASE (1)
       AM%AM => AMAG1
       AM%DESC = 'AMAG1'
    CASE DEFAULT
       ! should never happen
       AM%AM => NULL()
       AM%DESC = 'UNKNOWN!'
       INFO = 0
    END SELECT
  END SUBROUTINE AMP_INIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE FUNCTION JSTEP_LEN(N, N_2)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, N_2
    INTEGER :: JSTEP_LEN

    IF (N_2 .EQ. 0) THEN
       JSTEP_LEN = MAX((N / 2), 0)
    ELSE
       JSTEP_LEN = MIN(MAX(N_2, 0), (N / 2))
    END IF
  END FUNCTION JSTEP_LEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE INIT_TRIU(N, P, Q, NN, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, NN
    INTEGER, INTENT(OUT) :: P(NN), Q(NN), INFO

    INTEGER :: IP, IQ, I

    IF (N .LT. 0) THEN
       INFO = -1
    ELSE IF (NN .LT. 0) THEN
       INFO = -4
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN
    IF (N .LE. 1) RETURN

    INFO = (N * (N - 1)) / 2

    I = 1
    DO IQ = 2, N
       DO IP = 1, IQ-1
          IF (I .GT. NN) RETURN
          P(I) = IP
          Q(I) = IQ
          I = I + 1
       END DO
    END DO
  END SUBROUTINE INIT_TRIU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE BUILD_JSTEP(N, A, LDA, J, NN, P, Q, AM, DZ, N_2, S, INFO)
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
  END SUBROUTINE BUILD_JSTEP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE JSTEP
