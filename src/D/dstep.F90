MODULE DSTEP
  USE OMP_LIB
  USE DTRANSF
  USE DTYPES
  USE JSTEP
  IMPLICIT NONE

CONTAINS

#ifdef MKL_DIRECT_CALL_SEQ_JIT
#include "mkl_direct_call.fi"
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE REAL(KIND=DWP) FUNCTION DMAG1(N, P, Q, A, LDA, J)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, P, Q, LDA, J(N)
    REAL(KIND=DWP), INTENT(IN) :: A(LDA,N)

    IF ((A(Q,P) .NE. D_ZERO) .OR. (A(P,Q) .NE. D_ZERO) .OR. (SIGN(D_ONE, A(P,P)) .EQ. D_MONE) .OR. &
         (SIGN(D_ONE, A(Q,Q)) .EQ. D_MONE) .OR. ((J(P) .EQ. J(Q)) .AND. (A(P,P) .LT. A(Q,Q)))) THEN
       DMAG1 = ABS(A(Q,P)) + ABS(A(P,Q))
    ELSE ! no transform
       DMAG1 = QUIET_NAN((P - 1) * N + (Q - 1))
    END IF
  END FUNCTION DMAG1

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

  SUBROUTINE DPROC_INIT(ID_MAG, ID_CMP, ID_TRU, R, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: ID_MAG, ID_CMP, ID_TRU
    TYPE(DPROC), INTENT(OUT) :: R
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
       R%MAG => DMAG1
    CASE (2)
       R%MAG => DMAG2
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

    IF (INT(OMP_GET_MAX_THREADS()) .GT. 1) THEN
       R%SRT => AW_SRT1
       R%NCP => AW_NCP1
    ELSE ! <= 1 (single-threaded)
       R%SRT => AW_SRT2
       R%NCP => AW_NCP2
    END IF
  END SUBROUTINE DPROC_INIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DSTEP_BUILD(S, N, A, LDA, J, NN, P, Q, R, DZ, N_2, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: S, N, LDA, J(N), NN, P(NN), Q(NN), N_2
    REAL(KIND=DWP), INTENT(IN) :: A(LDA,N)
    TYPE(DPROC), INTENT(IN) :: R
    TYPE(AW), INTENT(OUT), TARGET :: DZ(NN)
    INTEGER, INTENT(OUT) :: SL, STEP(N_2), INFO

    INTEGER :: IP, IQ, I, II, IT
#ifndef NDEBUG
    REAL(KIND=DWP) :: T
#endif

    IF (S .LT. 0) THEN
       INFO = -1
    ELSE IF (N .LT. 0) THEN
       INFO = -2
    ELSE IF (LDA .LT. N) THEN
       INFO = -4
    ELSE IF (NN .LT. 0) THEN
       INFO = -6
    ELSE IF (N_2 .LT. 0) THEN
       INFO = -11
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
       DZ(I)%P = IP
       DZ(I)%Q = IQ
       DZ(I)%B = IQ - IP
       DZ(I)%W = R%MAG(N, IP, IQ, A, LDA, J)
       IF (DZ(I)%W .EQ. DZ(I)%W) IT = IT + 1
    END DO
    !$OMP END PARALLEL DO

    IF (IT .EQ. 0) GOTO 1
#ifndef NDEBUG
    I = OPEN_LOG('DSTEP_BUILD', S)
#endif
    CALL R%SRT(NN, DZ, R%CMP, II)
    IF (II .LT. 0) THEN
       INFO = -10
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
    CALL R%NCP(NN, DZ, IT, SL, STEP, II)
    IF (II .LT. 0) THEN
       INFO = -12
       RETURN
    END IF
#ifndef NDEBUG
    IF (I .NE. -1) THEN
       T = II * DNS2S
       WRITE (I,'(A,F12.6,A)',ADVANCE='NO') 'NCP: ', T, ' s, '
    END IF
#endif

1   INFO = MAX(GET_THREAD_NS() - INFO, 1)
#ifndef NDEBUG
    IF (I .NE. -1) THEN
       T = INFO * DNS2S
       WRITE (I,'(A,F12.6,A)') 'BUILD: ', T, ' s'
       CALL AW_OUT(I, '', NN, DZ, SL, STEP, II)
       CLOSE(UNIT=I, IOSTAT=II)
    END IF
#endif
  END SUBROUTINE DSTEP_BUILD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DSTEP_TRANSF(S, N, U, LDU, A, LDA, Z, LDZ, J, NN, P, Q, R, DZ, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: S, N, LDU, LDA, LDZ, J(N), NN, P(NN), Q(NN), SL
    INTEGER, INTENT(INOUT) :: STEP(SL)
    REAL(KIND=DWP), INTENT(INOUT) :: U(LDU,N), A(LDA,N), Z(LDZ,N)
    TYPE(DPROC), INTENT(IN) :: R
    TYPE(AW), INTENT(IN) :: DZ(NN)
    INTEGER, INTENT(OUT) :: INFO

    INFO = SL
  END SUBROUTINE DSTEP_TRANSF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DSTEP_EXEC(S, N, U, LDU, A, LDA, Z, LDZ, J, NN, P, Q, R, DZ, N_2, STEP, SL, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: S, N, LDU, LDA, LDZ, J(N), NN, P(NN), Q(NN), N_2
    REAL(KIND=DWP), INTENT(INOUT) :: U(LDU,N), A(LDA,N), Z(LDZ,N)
    TYPE(DPROC), INTENT(IN) :: R
    TYPE(AW), INTENT(OUT), TARGET :: DZ(NN)
    INTEGER, INTENT(OUT) :: STEP(N_2), SL, INFO

    INTEGER :: IT, NL

    SL = 0
    INFO = 0
    IF (CtrlC .NE. 0) RETURN

    WRITE (ULOG,'(I10,A)',ADVANCE='NO') S, ','
    FLUSH(ULOG)
    CALL DSTEP_BUILD(S, N, A, LDA, J, NN, P, Q, R, DZ, N_2, SL, STEP, INFO)
    IF (INFO .LE. 0) THEN
       NL = INFO
    ELSE
       NL = SL
    END IF
    WRITE (ULOG,'(I11,A)',ADVANCE='NO') NL, ','
    IF (INFO .LE. 0) THEN
       WRITE (ULOG,'(F12.6,A,F12.6)') D_MZERO, ',', D_ZERO
    ELSE IF (SL .LE. 0) THEN
       WRITE (ULOG,'(F12.6,A,F12.6)') (INFO * DNS2S), ',', D_ZERO
    ELSE
       WRITE (ULOG,'(F12.6,A)',ADVANCE='NO') (INFO * DNS2S), ','
    END IF
    FLUSH(ULOG)

    IF (SL .LE. 0) RETURN
    IF (CtrlC .NE. 0) INFO = 0
    IF (INFO .LE. 0) THEN
       WRITE (ULOG,'(F12.6)') D_MZERO
       FLUSH(ULOG)
       RETURN
    END IF

    IT = GET_THREAD_NS()
    CALL DSTEP_TRANSF(S, N, U, LDU, A, LDA, Z, LDZ, J, NN, P, Q, R, DZ, SL, STEP, NL)
    IT = MAX(GET_THREAD_NS() - IT, 1)
    WRITE (ULOG,'(F12.6)') (IT * DNS2S)
    FLUSH(ULOG)

    SL = MIN(SL, NL)
    INFO = INFO + IT
  END SUBROUTINE DSTEP_EXEC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DSTEP_LOOP(N, U, LDU, A, LDA, Z, LDZ, J, NN, P, Q, R, DZ, N_2, STEP, INFO)
#ifdef ANIMATE
    USE VN_MTXVIS_F
#endif
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, LDU, LDA, LDZ, J(N), NN, P(NN), Q(NN), N_2
    REAL(KIND=DWP), INTENT(INOUT) :: U(LDU,N), A(LDA,N), Z(LDZ,N)
    TYPE(DPROC), INTENT(IN) :: R
    TYPE(AW), INTENT(OUT), TARGET :: DZ(NN)
    INTEGER, INTENT(OUT) :: STEP(N_2), INFO

    INTEGER :: S, SL
#ifdef ANIMATE
    CHARACTER(LEN=7), PARAMETER :: FNAME = 'djkAstp'
    INTEGER, PARAMETER :: ACT = IOR(IOR(VN_MTXVIS_OP_A, VN_MTXVIS_FN_Lg), VN_MTXVIS_FF_Bin) !VN_MTXVIS_FN_Id
    INTEGER, PARAMETER :: SX = 1, SY = 1
    TYPE(c_ptr) :: CTX

    INFO = VN_MTXVIS_START(CTX, FNAME, ACT, N, N, SX, SY, LEN_TRIM(FNAME))
    IF (INFO .NE. 0) THEN
       WRITE (ULOG,'(A,I11)') 'VN_MTXVIS_START:', INFO
       RETURN
    END IF
#else
    INFO = 0
#endif

    WRITE (ULOG,'(A)') '"STEP","BUILDs","TRANSFs"'
    FLUSH(ULOG)

    S = 0
    DO WHILE ((S .GE. 0) .AND. (CtrlC .EQ. 0))
#ifdef ANIMATE
       INFO = VN_MTXVIS_FRAME(CTX, A, N)
       IF (INFO .NE. 0) THEN
          WRITE (ULOG,'(A,I11)') 'VN_MTXVIS_FRAME:', INFO
          RETURN
       END IF
#endif
       CALL DSTEP_EXEC(S, N, U, LDU, A, LDA, Z, LDZ, J, NN, P, Q, R, DZ, N_2, STEP, SL, INFO)
       IF (INFO .LE. 0) EXIT
       IF (SL .LE. 0) EXIT
       S = S + 1
    END DO
    IF (INFO .GE. 0) INFO = S

#ifdef ANIMATE
    S = VN_MTXVIS_STOP(CTX)
    IF (S .NE. 0) WRITE (ULOG,'(A,I11)') 'VN_MTXVIS_STOP:', S
#endif
  END SUBROUTINE DSTEP_LOOP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE DSTEP
