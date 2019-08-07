MODULE ZSTEP
  USE ZTYPES
  USE JSTEP
  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE FUNCTION ZMAG1(N, P, Q, A, LDA, J)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, P, Q, LDA, J(N)
    COMPLEX(KIND=DWP), INTENT(IN) :: A(LDA,N)
    REAL(KIND=DWP) :: ZMAG1

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
       ZMAG1 = AAQP + AAPQ
    ELSE ! no transform
       ZMAG1 = QUIET_NAN((P - 1) * N + (Q - 1))
    END IF
  END FUNCTION ZMAG1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE ZPROC_INIT(ID_MAG, ID_CMP, ID_TRU, R, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: ID_MAG, ID_CMP, ID_TRU
    TYPE(ZPROC), INTENT(OUT) :: R
    INTEGER, INTENT(OUT) :: INFO

    IF (ID_MAG .LT. 0) THEN
       INFO = -1
    ELSE IF (ID_CMP .LT. 0) THEN
       INFO = -2
    ELSE IF (ID_TRU .LT. 0) THEN
       INFO = -3
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    IF (ID_MAG .EQ. 0) ID_MAG = 1
    SELECT CASE (ID_MAG)
    CASE (1)
       R%MAG => ZMAG1
    CASE DEFAULT
       R%MAG => NULL()
       INFO = -1
    END SELECT

    IF (ID_CMP .EQ. 0) ID_CMP = 1
    SELECT CASE (ID_CMP)
    CASE (1)
       R%CMP => AW_CMP1
    CASE DEFAULT
       R%CMP => NULL()
       INFO = -2
    END SELECT

    IF (ID_TRU .EQ. 0) ID_TRU = 1
    SELECT CASE (ID_TRU)
    CASE (1)
       R%PQI => PQI1
       R%TRU => TRU1
    CASE (2)
       R%PQI => PQI2
       R%TRU => TRU2
    CASE DEFAULT
       R%TRU => NULL()
       INFO = -3
    END SELECT
  END SUBROUTINE ZPROC_INIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE BUILD_ZSTEP(N, A, LDA, J, NN, P, Q, R, DZ, N_2, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, LDA, J(N), NN, P(NN), Q(NN), N_2
    COMPLEX(KIND=DWP), INTENT(IN) :: A(LDA,N)
    TYPE(ZPROC), INTENT(IN) :: R
    TYPE(AW), INTENT(OUT), TARGET :: DZ(NN)
    INTEGER, INTENT(OUT) :: SL, STEP(N_2), INFO

    INTEGER :: IP, IQ, I, II, IT
#ifndef NDEBUG
    REAL(KIND=DWP) :: T
#endif

    IF (N .LT. 0) THEN
       INFO = -1
    ELSE IF (LDA .LT. N) THEN
       INFO = -3
    ELSE IF (NN .LT. 0) THEN
       INFO = -5
    ELSE IF (N_2 .LT. 0) THEN
       INFO = -10
    ELSE ! all OK
       !DIR$ VECTOR ALWAYS
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    INFO = GET_THREAD_NS()
    SL = 0
    IF (N .EQ. 0) GOTO 1
    IF (NN .EQ. 0) GOTO 1
    IF (N_2 .EQ. 0) GOTO 1

    IT = 0
    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(IP,IQ,I) SHARED(NN,N,A,LDA,J,P,Q,R,DZ) REDUCTION(+:IT)
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
       DZ(I)%W = R%MAG(N, IP, IQ, A, LDA, J)
       IF (DZ(I)%W .EQ. DZ(I)%W) IT = IT + 1
    END DO
    !$OMP END PARALLEL DO

    IF (IT .EQ. 0) GOTO 1
#ifndef NDEBUG
    I = OPEN_LOG('BUILD_JSTEP')
#endif
    CALL AW_SORT(NN, DZ, R%CMP, II)
    IF (II .LT. 0) THEN
       INFO = -9
       RETURN
    END IF
#ifndef NDEBUG
    IF (I .NE. -1) THEN
       T = II * DNS2S
       CALL AW_OUT(I, '', NN, DZ, 0, STEP, II)
       WRITE (I,'(A,F12.6,A)',ADVANCE='NO') 'SORT: ', T, ' s, '
    END IF
#endif

    IT = MIN(IT, N_2)
    CALL AW_NCP(NN, DZ, IT, SL, STEP, II)
    IF (II .LT. 0) THEN
       INFO = -11
       RETURN
    END IF
#ifndef NDEBUG
    IF (I .NE. -1) THEN
       T = II * DNS2S
       WRITE (I,'(A,F12.6,A)',ADVANCE='NO') 'NCP: ', T, ' s, '
    END IF
#endif

1   INFO = GET_THREAD_NS() - INFO
#ifndef NDEBUG
    IF (I .NE. -1) THEN
       T = INFO * DNS2S
       WRITE (I,'(A,F12.6,A)') 'BUILD: ', T, ' s'
       CALL AW_OUT(I, '', NN, DZ, SL, STEP, II)
       CLOSE(UNIT=I, IOSTAT=II)
    END IF
#endif
  END SUBROUTINE BUILD_ZSTEP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE ZSTEP
