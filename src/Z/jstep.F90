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

  PURE FUNCTION MAG1(APP, AQP, APQ, AQQ, JP, JQ)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(IN) :: APP, AQP, APQ, AQQ
    INTEGER, INTENT(IN) :: JP, JQ
    REAL(KIND=DWP) :: MAG1

    MAG1 = ABS(AQP) + ABS(APQ)
  END FUNCTION MAG1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE FUNCTION CVG1(APP, AQP, APQ, AQQ, JP, JQ)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(IN) :: APP, AQP, APQ, AQQ
    INTEGER, INTENT(IN) :: JP, JQ
    INTEGER :: CVG1

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
       CVG1 = 1
    ELSE ! no transform
       CVG1 = 0
    END IF
  END FUNCTION CVG1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! column-cyclic
  PURE SUBROUTINE TRU1(N, P, Q, NN, INFO)
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
  END SUBROUTINE TRU1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! row-cyclic
  PURE SUBROUTINE TRU2(N, P, Q, NN, INFO)
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
    DO IP = 1, N-1
       DO IQ = IP+1, N
          IF (I .GT. NN) RETURN
          P(I) = IP
          Q(I) = IQ
          I = I + 1
       END DO
    END DO
  END SUBROUTINE TRU2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE APROC_INIT(ID_MAG, ID_CMP, ID_CVG, ID_TRU, R, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: ID_MAG, ID_CMP, ID_CVG, ID_TRU
    TYPE(APROC), INTENT(OUT) :: R
    INTEGER, INTENT(OUT) :: INFO

    IF (ID_MAG .LT. 0) THEN
       INFO = -1
    ELSE IF (ID_CMP .LT. 0) THEN
       INFO = -2
    ELSE IF (ID_CVG .LT. 0) THEN
       INFO = -3
    ELSE IF (ID_TRU .LT. 0) THEN
       INFO = -4
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    IF (ID_MAG .EQ. 0) ID_MAG = 1
    SELECT CASE (ID_MAG)
    CASE (1)
       R%MAG => MAG1
    CASE DEFAULT
       ! should never happen
       R%MAG => NULL()
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
       R%CVG => CVG1
    CASE DEFAULT
       ! should never happen
       R%CVG => NULL()
       INFO = -3
    END SELECT

    IF (ID_TRU .EQ. 0) ID_TRU = 1
    SELECT CASE (ID_TRU)
    CASE (1)
       R%TRU => TRU1
    CASE (2)
       R%TRU => TRU2
    CASE DEFAULT
       ! should never happen
       R%TRU => NULL()
       INFO = -4
    END SELECT
  END SUBROUTINE APROC_INIT

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

  SUBROUTINE BUILD_JSTEP(N, A, LDA, J, NN, P, Q, R, DZ, N_2, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, LDA, J(N), NN, P(NN), Q(NN), N_2
    COMPLEX(KIND=DWP), INTENT(IN) :: A(LDA,N)
    TYPE(APROC), INTENT(IN) :: R
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
       IT = R%CVG(A(IP,IP), A(IQ,IP), A(IP,IQ), A(IQ,IQ), J(IP), J(IQ))
       IF (IT .NE. 0) THEN
          DZ(I)%W = R%MAG(A(IP,IP), A(IQ,IP), A(IP,IQ), A(IQ,IQ), J(IP), J(IQ))
       ELSE ! NaN
          DZ(I)%W = QUIET_NAN(I)
       END IF
    END DO
    !$OMP END PARALLEL DO

    IF (IT .EQ. 0) RETURN

    CALL DZBW_SORT(NN, DZ, R%CMP, INFO)
    IF (INFO .LT. 0) RETURN

    IT = MIN(IT, N_2)
    CALL DZBW_NCP(NN, DZ, IT, STEP, INFO)
  END SUBROUTINE BUILD_JSTEP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE JSTEP
