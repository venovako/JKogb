MODULE ZSTEP
  USE, INTRINSIC :: ISO_C_BINDING
  USE OMP_LIB
  USE ZTRANSF
  USE ZTYPES
  USE JSTEP
  USE UTILS
  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE REAL(KIND=DWP) FUNCTION ZMAG1(N, P, Q, A, LDA, J)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, P, Q, LDA, J(N)
    COMPLEX(KIND=DWP), INTENT(IN) :: A(LDA,N)

    IF ((A(Q,P) .NE. Z_ZERO) .OR. (A(P,Q) .NE. Z_ZERO) .OR. (AIMAG(A(P,P)) .NE. D_ZERO) .OR. (AIMAG(A(Q,Q)) .NE. D_ZERO) .OR. &
         (SIGN(D_ONE, REAL(A(P,P))) .EQ. D_MONE) .OR. (SIGN(D_ONE, REAL(A(Q,Q))) .EQ. D_MONE) .OR. &
         ((J(P) .EQ. J(Q)) .AND. (REAL(A(P,P)) .LT. REAL(A(Q,Q))))) THEN
       ZMAG1 = ABS(A(Q,P)) + ABS(A(P,Q))
    ELSE ! no transform
       ZMAG1 = QUIET_NAN((P - 1) * N + (Q - 1))
    END IF
  END FUNCTION ZMAG1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE ZPROC_INIT(NT, ID_MAG, ID_CMP, ID_TRU, R, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT
    INTEGER, INTENT(INOUT) :: ID_MAG, ID_CMP, ID_TRU
    TYPE(ZPROC), INTENT(OUT) :: R
    INTEGER, INTENT(OUT) :: INFO

    IF (NT .LE. 0) THEN
       INFO = -1
    ELSE IF (ID_MAG .LT. 0) THEN
       INFO = -2
    ELSE IF (ID_CMP .LT. 0) THEN
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
       R%MAG => ZMAG1
    CASE DEFAULT
       R%MAG => NULL()
       INFO = -2
    END SELECT

    IF (ID_CMP .EQ. 0) ID_CMP = 1
    SELECT CASE (ID_CMP)
    CASE (1)
       R%CMP => AW_CMP1
    CASE DEFAULT
       R%CMP => NULL()
       INFO = -3
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
       INFO = -4
    END SELECT

    IF (NT .GT. 1) THEN
       R%SRT => AW_SRT1
       R%NCP => AW_NCP1
    ELSE ! <= 1 (single-threaded)
       R%SRT => AW_SRT2
       R%NCP => AW_NCP2
    END IF
  END SUBROUTINE ZPROC_INIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE ZSTEP_BUILD(NT, S, N, A, LDA, J, NN, P, Q, R, NM, DZ, N_2, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, S, N, LDA, J(N), NN, P(NN), Q(NN), NM, N_2
    COMPLEX(KIND=DWP), INTENT(IN) :: A(LDA,N)
    TYPE(ZPROC), INTENT(IN) :: R
    TYPE(AW), INTENT(OUT), TARGET :: DZ(NM)
    INTEGER, INTENT(OUT) :: SL, STEP(N_2), INFO

    INTEGER :: IP, IQ, I, II, IT
#ifndef NDEBUG
    REAL(KIND=DWP) :: T
#endif

    IF (NT .LE. 0) THEN
       INFO = -1
    ELSE IF (S .LT. 0) THEN
       INFO = -2
    ELSE IF (N .LT. 0) THEN
       INFO = -3
    ELSE IF (LDA .LT. N) THEN
       INFO = -5
    ELSE IF (NN .LT. 0) THEN
       INFO = -7
    ELSE IF (NM .LT. NN) THEN
       INFO = -11
    ELSE IF (N_2 .LT. 0) THEN
       INFO = -13
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    INFO = GET_THREAD_NS()
    SL = 0
    IF (N .EQ. 0) GOTO 1
    IF (NN .EQ. 0) GOTO 1
    IF (N_2 .EQ. 0) GOTO 1

    IT = 0
    !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(IP,IQ,I) SHARED(NT,NN,N,A,LDA,J,P,Q,R,DZ) REDUCTION(+:IT)
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
    I = OPEN_LOG('ZSTEP_BUILD', S)
#endif
    CALL R%SRT(NT, NN, NM, DZ, R%CMP, II)
    IF (II .LT. 0) THEN
       INFO = -11
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
    CALL R%NCP(NT, NN, NM, DZ, IT, SL, STEP, II)
    IF (II .LT. 0) THEN
       INFO = -13
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
  END SUBROUTINE ZSTEP_BUILD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE ZSTEP_TRANSF(NT, S, N, U, LDU, A, LDA, Z, LDZ, J, SIGMA, NN, NM, DZ, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, S, N, LDU, LDA, LDZ, J(N), NN, NM, SL
    INTEGER, INTENT(INOUT) :: STEP(SL)
    COMPLEX(KIND=DWP), INTENT(INOUT) :: U(LDU,N), A(LDA,N), Z(LDZ,N)
    REAL(KIND=DWP), INTENT(OUT), TARGET :: SIGMA(N)
    TYPE(AW), INTENT(INOUT), TARGET :: DZ(NM)
    INTEGER, INTENT(OUT) :: INFO

    COMPLEX(KIND=DWP) :: V(2,2), A2(2,2)
    COMPLEX(KIND=DWP), POINTER, CONTIGUOUS :: W(:,:,:)
    INTEGER :: I, P, Q
    INTEGER, POINTER, CONTIGUOUS :: IT(:)

    CALL C_F_POINTER(C_LOC(DZ(NN+1)), W, [2,2,SL])
    CALL C_F_POINTER(C_LOC(SIGMA), IT, [SL])
    INFO = HUGE(0)

!$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(I,P,Q,V,A2) SHARED(N,SL,STEP,IT,DZ,J,U,A,Z,LDA,W) REDUCTION(MIN:INFO)
    DO I = 1, SL
       P = DZ(STEP(I))%P
       Q = DZ(STEP(I))%Q

       A2(1,1) = A(P,P)
       A2(2,1) = A(Q,P)
       A2(1,2) = A(P,Q)
       A2(2,2) = A(Q,Q)

       CALL ZHSVD2((J(P) .NE. J(Q)), A2, V, W(:,:,I), IT(I))
       INFO = IT(I)

       IF (IT(I) .GT. 1) THEN
          IF (IAND(IT(I), 2) .NE. 0) THEN
             CALL BA(V, N, A(P,1), A(Q,1), LDA)
             CALL UH(V)
             CALL AB(V, N, U(1,P), U(1,Q))
          END IF
          IF (IAND(IT(I), 4) .NE. 0) CALL AB(W(:,:,I), N, Z(1,P), Z(1,Q))
       END IF
    END DO
!$OMP END PARALLEL DO

    IF (INFO .GE. 0) THEN
       INFO = 0
       !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(I,P,Q) SHARED(N,SL,STEP,IT,DZ,A,LDA,W) REDUCTION(+:INFO)
       DO I = 1, SL
          P = DZ(STEP(I))%P
          Q = DZ(STEP(I))%Q

          IF (IAND(IT(I), 4) .NE. 0) CALL AB(W(:,:,I), N, A(1,P), A(1,Q))
          IF ((IAND(IT(I), 2) .NE. 0) .OR. (IAND(IT(I), 4) .NE. 0)) INFO = INFO + 1

          A(Q,P) = Z_ZERO
          A(P,Q) = Z_ZERO
       END DO
       !$OMP END PARALLEL DO
    END IF
  END SUBROUTINE ZSTEP_TRANSF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE ZSTEP_EXEC(NT, S, N, U, LDU, A, LDA, Z, LDZ, J, SIGMA, NN, P, Q, R, NM, DZ, N_2, STEP, SL, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, S, N, LDU, LDA, LDZ, J(N), NN, P(NN), Q(NN), NM, N_2
    COMPLEX(KIND=DWP), INTENT(INOUT) :: U(LDU,N), A(LDA,N), Z(LDZ,N)
    REAL(KIND=DWP), INTENT(OUT), TARGET :: SIGMA(N)
    TYPE(ZPROC), INTENT(IN) :: R
    TYPE(AW), INTENT(OUT), TARGET :: DZ(NM)
    INTEGER, INTENT(OUT) :: STEP(N_2), SL, INFO

    INTEGER :: IT, NL

    SL = 0
    INFO = 0
    IF (CtrlC .NE. 0_ATOMIC_INT_KIND) RETURN

    WRITE (ULOG,'(I10,A)',ADVANCE='NO') S, ','
    FLUSH(ULOG)
    CALL ZSTEP_BUILD(NT, S, N, A, LDA, J, NN, P, Q, R, NM, DZ, N_2, SL, STEP, INFO)
    IF (INFO .LE. 0) THEN
       IT = INFO
    ELSE
       IT = SL
    END IF
    WRITE (ULOG,'(I11,A)',ADVANCE='NO') IT, ','
    FLUSH(ULOG)

    IF (INFO .LE. 0) THEN
       WRITE (ULOG,'(F12.6,A,F12.6,A,I11)') D_MZERO, ',', D_ZERO, ',', -1
    ELSE IF (SL .LE. 0) THEN
       WRITE (ULOG,'(F12.6,A,F12.6,A,I11)') (INFO * DNS2S), ',', D_ZERO, ',', -1
    ELSE
       WRITE (ULOG,'(F12.6,A)',ADVANCE='NO') (INFO * DNS2S), ','
    END IF
    FLUSH(ULOG)

    IF (SL .LE. 0) RETURN
    IF (CtrlC .NE. 0_ATOMIC_INT_KIND) INFO = 0
    IF (INFO .LE. 0) THEN
       WRITE (ULOG,'(F12.6,A,I11)') D_MZERO, ',', -1
       FLUSH(ULOG)
       RETURN
    END IF

    IT = GET_THREAD_NS()
    CALL ZSTEP_TRANSF(NT, S, N, U, LDU, A, LDA, Z, LDZ, J, SIGMA, NN, NM, DZ, SL, STEP, NL)
    IT = MAX(GET_THREAD_NS() - IT, 1)
    WRITE (ULOG,'(F12.6,A,I11)') (IT * DNS2S), ',', NL
    FLUSH(ULOG)

    SL = MIN(SL, NL)
    INFO = INFO + IT
  END SUBROUTINE ZSTEP_EXEC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE JZHJ(NT, N, Z, LDZ, J, NN, P, Q)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, N, LDZ, J(N), NN, P(NN), Q(NN)
    COMPLEX(KIND=DWP), INTENT(INOUT) :: Z(LDZ,N)

    COMPLEX(KIND=DWP) :: T
    INTEGER :: IP, IQ, I

    !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(IP,IQ,I,T) SHARED(NN,Z,J,P,Q)
    DO I = 1, NN
       IP = P(I)
       IQ = Q(I)
       T = Z(IP,IQ)
       Z(IP,IQ) = CONJG(Z(IQ,IP))
       Z(IQ,IP) = CONJG(T)
       IF (J(IP) .NE. J(IQ)) THEN
          Z(IP,IQ) = -Z(IP,IQ)
          Z(IQ,IP) = -Z(IQ,IP)
       END IF
    END DO
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I) SHARED(N,Z)
    DO I = 1, N
       Z(I,I) = CONJG(Z(I,I))
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE JZHJ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE ZSTEP_LOOP(NT, N, U, LDU, A, LDA, Z, LDZ, J, SIGMA, NN, P, Q, R, NM, DZ, N_2, STEP, INFO)
#ifdef ANIMATE
    USE VN_CMPLXVIS_F
#endif
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, N, LDU, LDA, LDZ, J(N), NN, P(NN), Q(NN), NM, N_2
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(LDA,N)
    COMPLEX(KIND=DWP), INTENT(OUT) :: U(LDU,N), Z(LDZ,N)
    REAL(KIND=DWP), INTENT(OUT), TARGET :: SIGMA(N)
    TYPE(ZPROC), INTENT(IN) :: R
    TYPE(AW), INTENT(OUT), TARGET :: DZ(NM)
    INTEGER, INTENT(OUT) :: STEP(N_2), INFO

    INTEGER :: S, SL
#ifdef ANIMATE
    CHARACTER(LEN=7), PARAMETER :: FNAME = 'zjkAstp'
    INTEGER, PARAMETER :: ACT = IOR(VN_CMPLXVIS_OP_A, VN_CMPLXVIS_FN_Lg) !VN_CMPLXVIS_FN_Id
    INTEGER, PARAMETER :: SX = ANIMATE, SY = ANIMATE
    TYPE(c_ptr) :: CTX
#endif

    IF (NT .LE. 0) THEN
       INFO = -1
    ELSE IF (N .LT. 3) THEN
       INFO = -2
    ELSE IF (LDU .LT. N) THEN
       INFO = -4
    ELSE IF (LDA .LT. N) THEN
       INFO = -6
    ELSE IF (LDZ .LT. N) THEN
       INFO = -8
    ELSE IF (NN .LT. 0) THEN
       INFO = -11
    ELSE IF (NM .LT. NN) THEN
       INFO = -15
    ELSE IF (N_2 .LT. 0) THEN
       INFO = -17
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(S,INFO) SHARED(N,U)
    DO S = 1, N
       !DIR$ VECTOR ALWAYS ASSERT
       DO INFO = 1, S-1
          U(INFO,S) = Z_ZERO
       END DO
       U(S,S) = Z_ONE
       !DIR$ VECTOR ALWAYS ASSERT
       DO INFO = S+1, N
          U(INFO,S) = Z_ZERO
       END DO
    END DO
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(S,INFO) SHARED(N,Z)
    DO S = 1, N
       !DIR$ VECTOR ALWAYS ASSERT
       DO INFO = 1, S-1
          Z(INFO,S) = Z_ZERO
       END DO
       Z(S,S) = Z_ONE
       !DIR$ VECTOR ALWAYS ASSERT
       DO INFO = S+1, N
          Z(INFO,S) = Z_ZERO
       END DO
    END DO
    !$OMP END PARALLEL DO

#ifdef ANIMATE
    INFO = VN_CMPLXVIS_START(CTX, FNAME, ACT, N, N, SX, SY, LEN_TRIM(FNAME))
    IF (INFO .NE. 0) THEN
       WRITE (ULOG,'(A,I11)') 'VN_CMPLXVIS_START:', INFO
       RETURN
    END IF
#else
    INFO = 0
#endif

    WRITE (ULOG,'(A)') '"STEP","OLDLEN","BUILDs","TRANSFs","NEWLEN"'
    FLUSH(ULOG)

    S = 0
    DO WHILE ((S .GE. 0) .AND. (CtrlC .EQ. 0_ATOMIC_INT_KIND))
#ifdef ANIMATE
       INFO = VN_CMPLXVIS_FRAME(CTX, A, N)
       IF (INFO .NE. 0) THEN
          WRITE (ULOG,'(A,I11)') 'VN_CMPLXVIS_FRAME:', INFO
          RETURN
       END IF
#endif
       CALL ZSTEP_EXEC(NT, S, N, U, LDU, A, LDA, Z, LDZ, J, SIGMA, NN, P, Q, R, NM, DZ, N_2, STEP, SL, INFO)
       IF (INFO .LE. 0) EXIT
       IF (SL .LE. 0) EXIT
       S = S + 1
    END DO
    IF (SL .LT. 0) THEN
       INFO = SL
    ELSE IF (INFO .GE. 0) THEN
       INFO = S
    END IF

#ifdef ANIMATE
    S = VN_CMPLXVIS_STOP(CTX)
    IF (S .NE. 0) WRITE (ULOG,'(A,I11)') 'VN_CMPLXVIS_STOP:', S
#endif

    !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(S) SHARED(N,SIGMA,A)
    DO S = 1, N
       SIGMA(S) = REAL(A(S,S))
       ! for a simple calculation of ||off(A)||_F
       A(S,S) = Z_ZERO
    END DO
    !$OMP END PARALLEL DO
    CALL JZHJ(NT, N, Z, LDZ, J, NN, P, Q)
  END SUBROUTINE ZSTEP_LOOP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE ZSTEP
