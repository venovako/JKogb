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

  PURE FUNCTION ACVG1(N, APP, AQP, APQ, AQQ, JP, JQ)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(IN) :: APP, AQP, APQ, AQQ
    INTEGER, INTENT(IN) :: N, JP, JQ
    INTEGER :: ACVG1

    REAL(KIND=DWP) :: AAPP, AAQQ, MAXPQ, MINPQ

    AAPP = ABS(APP)
    AAQQ = ABS(AQQ)
    IF (AAPP .GE. AAQQ) THEN
       MAXPQ = AAPP
       MINPQ = AAQQ
    ELSE
       MAXPQ = AAQQ
       MINPQ = AAPP
    END IF

    IF (MAX(ABS(AQP), ABS(APQ)) .GT. ((MAXPQ * D_EPS) * MINPQ)) THEN
       ACVG1 = 1
    ELSE ! no transform
       ACVG1 = 0
    END IF
  END FUNCTION ACVG1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE AMC_INIT(ID_AMP, ID_CMP, ID_CVG, R, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: ID_AMP, ID_CMP, ID_CVG
    TYPE(AMC), INTENT(OUT) :: R
    INTEGER, INTENT(OUT) :: INFO

    IF (ID_AMP .LT. 0) THEN
       INFO = -1
    ELSE IF (ID_CMP .LT. 0) THEN
       INFO = -2
    ELSE IF (ID_CVG .LT. 0) THEN
       INFO = -3
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    IF (ID_AMP .EQ. 0) ID_AMP = 1
    SELECT CASE (ID_AMP)
    CASE (1)
       R%AMP => AMAG1
    CASE DEFAULT
       ! should never happen
       R%AMP => NULL()
       INFO = -1
    END SELECT

    IF (ID_CMP .EQ. 0) ID_CMP = 1
    SELECT CASE (ID_CMP)
    CASE (1)
       R%CMP => DZBW_CMP1
    CASE (2)
       R%CMP => DZBW_CMP2
    CASE DEFAULT
       ! should never happen
       R%CMP => NULL()
       INFO = -2
    END SELECT

    IF (ID_CVG .EQ. 0) ID_CVG = 1
    SELECT CASE (ID_CVG)
    CASE (1)
       R%CVG => ACVG1
    CASE DEFAULT
       ! should never happen
       R%CVG => NULL()
       INFO = -3
    END SELECT
  END SUBROUTINE AMC_INIT

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

  SUBROUTINE BUILD_JSTEP(N, A, LDA, J, NN, P, Q, R, DZ, N_2, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, LDA, J(N), NN, P(NN), Q(NN), N_2
    COMPLEX(KIND=DWP), INTENT(IN) :: A(LDA,N)
    TYPE(AMC), INTENT(IN) :: R
    TYPE(DZBW), INTENT(OUT), TARGET :: DZ(NN)
    INTEGER, INTENT(OUT) :: STEP(N_2), INFO

    INTEGER :: IP, IQ, I, IT

    IF (N .LT. 0) THEN
       INFO = -1
    ELSE IF (LDA .LT. N) THEN
       INFO = -3
    ELSE IF (NN .LT. 0) THEN
       INFO = -5
    ELSE IF (N_2 .LT. 0) THEN
       INFO = -10
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN
    IF (N .EQ. 0) RETURN
    IF (NN .EQ. 0) RETURN
    IF (N_2 .EQ. 0) RETURN

    IT = 0
    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(IP,IQ,I) SHARED(NN,N,A,J,P,Q,R,DZ) REDUCTION(+:IT)
    DO I = 1, NN
       IP = P(I)
       IQ = Q(I)
       IF (J(IP) .GE. 0) THEN
          IF (J(IQ) .GE. 0) THEN
             DZ(I)%P = IP
             DZ(I)%Q = IQ
             DZ(I)%B = IQ - IP
          ELSE ! JQ < 0
             DZ(I)%P = IP
             DZ(I)%Q = -IQ
             DZ(I)%B = IP - IQ
          END IF
       ELSE ! JP < 0
          IF (J(IQ) .GE. 0) THEN
             DZ(I)%P = -IP
             DZ(I)%Q = IQ
             DZ(I)%B = IP - IQ
          ELSE ! JQ < 0
             DZ(I)%P = -IP
             DZ(I)%Q = -IQ
             DZ(I)%B = IQ - IP
          END IF
       END IF
       IT = R%CVG(N, A(IP,IP), A(IQ,IP), A(IP,IQ), A(IQ,IQ), J(IP), J(IQ))
       IF (IT .NE. 0) THEN
          DZ(I)%W = R%AMP(A(IP,IP), A(IQ,IP), A(IP,IQ), A(IQ,IQ), J(IP), J(IQ))
       ELSE ! NaN
          DZ(I)%W = TRANSFER(-1_DWP, D_ZERO)
       END IF
    END DO
    !$OMP END PARALLEL DO

    IF (IT .EQ. 0) RETURN

    CALL DZBW_SORT(NN, DZ, R%CMP, INFO)
    IF (INFO .NE. 0) RETURN

    IT = MIN(IT, N_2)
    CALL DZBW_NCP(NN, DZ, IT, STEP, INFO)
  END SUBROUTINE BUILD_JSTEP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE JSTEP
