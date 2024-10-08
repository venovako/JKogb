MODULE zstep
  !$ USE omp_lib
  USE ztransf
  USE ztypes
  USE jstep
  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE REAL(KIND=DWP) FUNCTION ZMAGF2T(N, P, Q, A, LDA, J)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, P, Q, LDA, J(N)
    COMPLEX(KIND=DWP), INTENT(IN) :: A(LDA,N)

    SELECT CASE (J(P) + J(Q))
    CASE (-2)
       IF ((A(Q,P) .NE. Z_ZERO) .OR. (A(P,Q) .NE. Z_ZERO) .OR. &
            (AIMAG(A(P,P)) .NE. D_ZERO) .OR. (AIMAG(A(Q,Q)) .NE. D_ZERO) .OR. &
            (SIGN(D_ONE, REAL(A(P,P))) .EQ. D_MONE) .OR. (SIGN(D_ONE, REAL(A(Q,Q))) .EQ. D_MONE) .OR. &
            (REAL(A(Q,Q)) .LT. REAL(A(P,P)))) THEN
          ZMAGF2T = DASUM4(REAL(A(Q,P)), AIMAG(A(Q,P)), REAL(A(P,Q)), AIMAG(A(P,Q)))
       ELSE ! no transform
          ZMAGF2T = QUIET_NAN((P - 1) * N + (Q - 1))
       END IF
    CASE (2)
       IF ((A(Q,P) .NE. Z_ZERO) .OR. (A(P,Q) .NE. Z_ZERO) .OR. &
            (AIMAG(A(P,P)) .NE. D_ZERO) .OR. (AIMAG(A(Q,Q)) .NE. D_ZERO) .OR. &
            (SIGN(D_ONE, REAL(A(P,P))) .EQ. D_MONE) .OR. (SIGN(D_ONE, REAL(A(Q,Q))) .EQ. D_MONE) .OR. &
            (REAL(A(P,P)) .LT. REAL(A(Q,Q)))) THEN
          ZMAGF2T = DASUM4(REAL(A(Q,P)), AIMAG(A(Q,P)), REAL(A(P,Q)), AIMAG(A(P,Q)))
       ELSE ! no transform
          ZMAGF2T = QUIET_NAN((P - 1) * N + (Q - 1))
       END IF
    CASE DEFAULT ! invalid J
       ZMAGF2T = QUIET_NAN((P - 1) * N + (Q - 1))
    END SELECT
  END FUNCTION ZMAGF2T

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE REAL(KIND=DWP) FUNCTION ZMAGF2H(N, P, Q, A, LDA, J)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, P, Q, LDA, J(N)
    COMPLEX(KIND=DWP), INTENT(IN) :: A(LDA,N)

    COMPLEX(KIND=DWP) :: A2(2,2), U2(2,2), Z2(2,2)
    INTEGER :: J2(2), INFO

    IF ((A(Q,P) .NE. Z_ZERO) .OR. (A(P,Q) .NE. Z_ZERO) .OR. (AIMAG(A(P,P)) .NE. D_ZERO) .OR. (AIMAG(A(Q,Q)) .NE. D_ZERO) .OR. &
         (SIGN(D_ONE, REAL(A(P,P))) .EQ. D_MONE) .OR. (SIGN(D_ONE, REAL(A(Q,Q))) .EQ. D_MONE)) THEN
       A2(1,1) = A(P,P)
       A2(2,1) = A(Q,P)
       A2(1,2) = A(P,Q)
       A2(2,2) = A(Q,Q)
       J2(1) = J(P)
       J2(2) = J(Q)
       CALL ZHSVD2(A2, J2, U2, Z2, TH1FIX, INFO)
       IF (INFO .LT. 0) THEN
          ZMAGF2H = QUIET_NAN((P - 1) * N + (Q - 1))
       ELSE IF (IAND(INFO, 1) .EQ. 0) THEN
          ZMAGF2H = D_ZERO
       ELSE IF (IAND(INFO, 4) .EQ. 0) THEN
          ZMAGF2H = DASUM4(REAL(A(Q,P)), AIMAG(A(Q,P)), REAL(A(P,Q)), AIMAG(A(P,Q)))
       ELSE IF (IAND(INFO, 8) .EQ. 0) THEN
          ZMAGF2H = ABODNZF2(Z2, N, A(1,P), A(1,Q), Z_ZERO, Z_ZERO, P, Q)
       ELSE ! a non-trivial transform
          ZMAGF2H = ABODNZF2(Z2, N, A(1,P), A(1,Q), A2(2,1), A2(1,2), P, Q)
       END IF
    ELSE ! no transform
       ZMAGF2H = QUIET_NAN((P - 1) * N + (Q - 1))
    END IF
  END FUNCTION ZMAGF2H

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE ZPROC_INIT(NT, ID_NCP, ID_TRU, R, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT
    INTEGER, INTENT(INOUT) :: ID_NCP, ID_TRU
    TYPE(ZPROC), INTENT(OUT) :: R
    INTEGER, INTENT(OUT) :: INFO

#ifdef USE_FAST
    INFO = 0
#else
    IF (NT .LE. 0) THEN
       INFO = -1
    ELSE IF (ID_NCP .LT. 0) THEN
       INFO = -2
    ELSE IF (ID_TRU .LT. 0) THEN
       INFO = -3
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN
#endif

    IF (ID_NCP .EQ. 0) THEN
       IF (NT .GT. 1) THEN
          ID_NCP = 1
       ELSE ! <= 1 (single-threaded)
          ID_NCP = 2
       END IF
    END IF

    SELECT CASE (ID_NCP)
    CASE (1)
       R%NCP => AW_NCP1
    CASE (2)
       R%NCP => AW_NCP2
    CASE (3)
       R%NCP => AW_NCP3
    CASE DEFAULT
       R%NCP => NULL()
       INFO = -2
    END SELECT

    IF (ID_TRU .EQ. 0) ID_TRU = 1
    SELECT CASE (ID_TRU)
    CASE (1)
       R%TRU => TRU1
    CASE (2)
       R%TRU => TRU2
    CASE DEFAULT
       R%TRU => NULL()
       INFO = -3
    END SELECT

    IF (NT .GT. 1) THEN
       R%SRT => AW_SRT1
    ELSE ! <= 1 (single-threaded)
       R%SRT => AW_SRT2
    END IF
  END SUBROUTINE ZPROC_INIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE ZSTEP_BUILD(S, N, A, LDA, J, NN, P, Q, R, NM, DZ, TT, N_2, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: S, N, LDA, J(N), NN, P(NN), Q(NN), NM, TT, N_2
    COMPLEX(KIND=DWP), INTENT(IN) :: A(LDA,N)
    TYPE(ZPROC), INTENT(IN) :: R
    TYPE(AW), INTENT(OUT) :: DZ(NM)
    INTEGER, INTENT(OUT) :: SL, STEP(N_2)
    INTEGER(KIND=DWP), INTENT(OUT) :: INFO

    TYPE(AW) :: X
    INTEGER(KIND=DWP) :: ID
    INTEGER :: IP, IQ, I
    INTEGER ::II, IT
    !$ SAVE ::II, IT

    SL = 0
#ifndef NDEBUG
    I = -1
#endif

#ifdef USE_FAST
    INFO = 0_DWP
#else
    IF (S .LT. 0) THEN
       INFO = -1_DWP
    ELSE IF (N .LT. 0) THEN
       INFO = -2_DWP
    ELSE IF (LDA .LT. N) THEN
       INFO = -4_DWP
    ELSE IF (NN .LT. 0) THEN
       INFO = -6_DWP
    ELSE IF (NM .LT. NN) THEN
       INFO = -10_DWP
    ELSE IF (TT .LT. 0) THEN
       INFO = -12_DWP
    ELSE IF (TT .GT. NN) THEN
       INFO = -12_DWP
    ELSE IF (N_2 .LT. 0) THEN
       INFO = -13_DWP
    ELSE ! all OK
       INFO = 0_DWP
    END IF
    IF (INFO .NE. 0_DWP) RETURN
    INFO = GET_SYS_TIME()
#endif

    IF (N .EQ. 0) GOTO 1
    IF (NN .EQ. 0) GOTO 1
    IF (N_2 .EQ. 0) GOTO 1

    IT = 0
    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(IP,IQ,I) SHARED(TT,N,A,LDA,J,P,Q,R,DZ) REDUCTION(+:IT)
    DO I = 1, TT
       IP = P(I)
       IQ = Q(I)
       DZ(I)%P = IP
       DZ(I)%Q = IQ
       DZ(I)%B = IQ - IP
       DZ(I)%I = I
       DZ(I)%W = ZMAGF2T(N, IP, IQ, A, LDA, J)
       IF (DZ(I)%W .EQ. DZ(I)%W) IT = IT + 1
    END DO
    !$OMP END PARALLEL DO

    II = 0
    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(IP,IQ,I) SHARED(NN,TT,N,A,LDA,J,P,Q,R,DZ) REDUCTION(+:II)
    DO I = TT+1, NN
       IP = P(I)
       IQ = Q(I)
       DZ(I)%P = IP
       DZ(I)%Q = IQ
       DZ(I)%B = IQ - IP
       DZ(I)%I = I
       DZ(I)%W = ZMAGF2H(N, IP, IQ, A, LDA, J)
       IF (DZ(I)%W .LT. -HUGE(DZ(I)%W)) THEN
          ! -\infty removed
          DZ(I)%W = QUIET_NAN(I)
       ELSE IF (DZ(I)%W .EQ. DZ(I)%W) THEN
          II = II + 1
       END IF
    END DO
    !$OMP END PARALLEL DO

    IT = IT + II
    IF (IT .EQ. 0) GOTO 1
    IF (IT .LT. NN) THEN
       ! remove NaN weights
       I = 1
       II = NN
       DO WHILE (I .LT. II)
          DO WHILE (.NOT. (DZ(II)%W .EQ. DZ(II)%W))
             II = II - 1
             IF (II .LE. I) EXIT
          END DO
          DO WHILE (I .LT. II)
             IF (DZ(I)%W .EQ. DZ(I)%W) THEN
                I = I + 1
             ELSE ! NaN
                X = DZ(I)
                DZ(I) = DZ(II)
                DZ(II) = X
                II = II - 1
                I = I + 1
                EXIT
             END IF
          END DO
       END DO
    END IF
#ifndef NDEBUG
    I = OPEN_LOG('ZSTEP_BUILD', S)
#endif
    CALL R%SRT(IT, NM, DZ, ID)
    IF (ID .LT. 0_DWP) THEN
       INFO = -11_DWP
       RETURN
    END IF
#ifndef NDEBUG
    IF (I .NE. -1) THEN
       WRITE (I,'(A,F12.6,A)') 'SRT: ', (ID / REAL(GET_SYS_TRES(),DWP)), ' s'
       CALL AW_OUT(I, '', IT, DZ, 0, STEP, ID)
    END IF
#endif

    CALL R%NCP(N, IT, NM, DZ, MIN(IT, N_2), SL, STEP, ID)
    IF (ID .LT. 0_DWP) THEN
       INFO = -13_DWP
       RETURN
    END IF
#ifndef NDEBUG
    IF (I .NE. -1) WRITE (I,'(A,F12.6,A)') 'NCP: ', (ID / REAL(GET_SYS_TRES(),DWP)), ' s'
#endif

1   CONTINUE
#ifndef USE_FAST
    INFO = GET_SYS_TLAP(INFO)
#endif

#ifndef NDEBUG
    IF (I .NE. -1) THEN
       WRITE (I,'(A,F12.6,A)') 'ALL: ', (INFO / REAL(GET_SYS_TRES(),DWP)), ' s'
       CALL AW_OUT(I, '', IT, DZ, SL, STEP, ID)
       CLOSE(UNIT=I, IOSTAT=II)
    END IF
#endif
  END SUBROUTINE ZSTEP_BUILD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE ZSTEP_TRANSF(N, U, LDU, A, LDA, Z, LDZ, J, SIGMA, NN, DZ, SL, STEP, W, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, LDU, LDA, LDZ, J(N), NN, SL, STEP(SL)
    COMPLEX(KIND=DWP), INTENT(INOUT) :: U(LDU,N), A(LDA,N), Z(LDZ,N)
    REAL(KIND=DWP), INTENT(OUT) :: SIGMA(N)
    TYPE(AW), INTENT(INOUT) :: DZ(NN)
    INTEGER, INTENT(OUT) :: INFO
    COMPLEX(KIND=DWP), INTENT(OUT) :: W(2,3,SL)

    COMPLEX(KIND=DWP) :: V(2,2), B(2,2)
    INTEGER :: I, L, K(2)
    REAL(KIND=DWP) :: TW
    INTEGER :: M, P, Q
    !$ SAVE :: TW, M, P, Q

    M = HUGE(M)
    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,P,Q,V,B,K) SHARED(N,SL,STEP,DZ,J,U,A,Z,LDA,W) REDUCTION(MIN:M)
    DO I = 1, SL
       P = DZ(STEP(I))%P
       Q = DZ(STEP(I))%Q

       B(1,1) = A(P,P)
       B(2,1) = A(Q,P)
       B(1,2) = A(P,Q)
       B(2,2) = A(Q,Q)
       K(1) = J(P)
       K(2) = J(Q)

       CALL ZHSVD2(B, K, V, W(:,:,I), TH1FIX, M)
       DZ(STEP(I))%I = M
       IF (M .GE. 0) THEN
          W(1,3,I) = REAL(B(1,1))
          W(2,3,I) = REAL(B(2,2))
          IF (IAND(M, 2) .NE. 0) THEN
             CALL BA(V, N, A(P,1), A(Q,1), LDA)
             CALL UH2(V)
             CALL AB(V, N, U(1,P), U(1,Q))
          END IF
          IF (IAND(M, 4) .NE. 0) CALL AB(W(:,:,I), N, Z(1,P), Z(1,Q))
       END IF
    END DO
    !$OMP END PARALLEL DO

    IF (M .GE. 0) THEN
       M = 0
       !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,L,P,Q) SHARED(N,SL,STEP,DZ,A,LDA,W) REDUCTION(+:M)
       DO I = 1, SL
          P = DZ(STEP(I))%P
          Q = DZ(STEP(I))%Q
          L = DZ(STEP(I))%I

          IF (IAND(L, 4) .NE. 0) CALL AB(W(:,:,I), N, A(1,P), A(1,Q))
          IF ((IAND(L, 1) .NE. 0) .AND. ((IAND(L, 2) .NE. 0) .OR. (IAND(L, 4) .NE. 0))) M = M + 1

          IF (IAND(L, 8) .EQ. 0) THEN
             A(P,P) = W(1,3,I)
             A(Q,P) = Z_ZERO
             A(P,Q) = Z_ZERO
             A(Q,Q) = W(2,3,I)
          END IF
       END DO
       !$OMP END PARALLEL DO
    END IF

    TW = D_MZERO
    P = 0
    Q = 0
    IF (M .GE. 0) THEN
       !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,L) SHARED(SL,STEP,DZ,J) REDUCTION(+:TW,P,Q)
       DO I = 1, SL
          L = DZ(STEP(I))%I
          IF ((IAND(L, 2) .NE. 0) .OR. (IAND(L, 4) .NE. 0)) THEN
             TW = TW + DZ(STEP(I))%W
             IF (J(DZ(STEP(I))%P) .EQ. J(DZ(STEP(I))%Q)) THEN
                P = P + 1
             ELSE ! hyperbolic
                Q = Q + 1
             END IF
          END IF
       END DO
       !$OMP END PARALLEL DO
    END IF

    SIGMA(1) = TW
    SIGMA(2) = P
    SIGMA(3) = Q
    INFO = M
  END SUBROUTINE ZSTEP_TRANSF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE ZSTEP_EXEC(S, N, U, LDU, A, LDA, Z, LDZ, J, SIGMA, NN, P, Q, R, NM, DZ, TT, N_2, STEP, SL, W, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: S, N, LDU, LDA, LDZ, J(N), NN, P(NN), Q(NN), NM, TT, N_2
    COMPLEX(KIND=DWP), INTENT(INOUT) :: U(LDU,N), A(LDA,N), Z(LDZ,N)
    REAL(KIND=DWP), INTENT(OUT) :: SIGMA(N)
    TYPE(ZPROC), INTENT(IN) :: R
    TYPE(AW), INTENT(OUT) :: DZ(NM)
    INTEGER, INTENT(OUT) :: STEP(N_2), SL, INFO
    COMPLEX(KIND=DWP), INTENT(OUT) :: W(2,3,N_2)

    INTEGER(KIND=DWP) :: IB, IT
    INTEGER :: NL

#ifndef USE_FAST
    WRITE (ERROR_UNIT,'(I11,A)',ADVANCE='NO') S, ','
    FLUSH(ERROR_UNIT)
#endif

    SL = 0
    CALL ZSTEP_BUILD(S, N, A, LDA, J, NN, P, Q, R, NM, DZ, TT, N_2, SL, STEP, IB)
    IF (IB .LT. 0_DWP) THEN
       INFO = INT(IB)
    ELSE ! IB .GE. 0
       INFO = 0
    END IF

#ifndef USE_FAST
    IF (INFO .LT. 0) THEN
       WRITE (ERROR_UNIT,'(I11,A)',ADVANCE='NO') INFO, ','
    ELSE ! INFO .GE. 0
       WRITE (ERROR_UNIT,'(I11,A)',ADVANCE='NO') SL, ','
    END IF
    FLUSH(ERROR_UNIT)

    IF (INFO .LT. 0) THEN
       WRITE (ERROR_UNIT,'(F12.6,A,F12.6,A,I11,A,ES25.17E3,2(A,I11))') &
            D_MZERO, ',', D_ZERO, ',', -1, ',', D_MZERO, ',',0,',',0
    ELSE IF (SL .LE. 0) THEN
       WRITE (ERROR_UNIT,'(F12.6,A,F12.6,A,I11,A,ES25.17E3,2(A,I11))') &
            (IB / REAL(GET_SYS_TRES(),DWP)), ',', D_ZERO, ',', -1, ',', D_MZERO, ',',0,',',0
    ELSE ! default
       WRITE (ERROR_UNIT,'(F12.6,A)',ADVANCE='NO') (IB / REAL(GET_SYS_TRES(),DWP)), ','
    END IF
    FLUSH(ERROR_UNIT)
#endif

    IF (SL .LE. 0) RETURN
    IF (INFO .LT. 0) THEN
#ifndef USE_FAST
       WRITE (ERROR_UNIT,'(F12.6,A,I11,A,ES25.17E3,2(A,I11))') &
            D_MZERO, ',', -1, ',', D_MZERO, ',',0,',',0
       FLUSH(ERROR_UNIT)
#endif
       RETURN
    END IF

#ifndef USE_FAST
    IT = GET_SYS_TIME()
#endif
    CALL ZSTEP_TRANSF(N, U, LDU, A, LDA, Z, LDZ, J, SIGMA, NN, DZ, SL, STEP, W, NL)
#ifndef USE_FAST
    IT = GET_SYS_TLAP(IT)
    WRITE (ERROR_UNIT,'(F12.6,A,I11,A,ES25.17E3,2(A,I11))') &
         (IT / REAL(GET_SYS_TRES(),DWP)), ',', NL, ',', SIGMA(1), ',',INT(SIGMA(2)),',',INT(SIGMA(3))
    FLUSH(ERROR_UNIT)
#endif

    SL = MIN(SL, NL)
  END SUBROUTINE ZSTEP_EXEC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE JZHJ(N, Z, LDZ, J, NN, P, Q)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, LDZ, J(N), NN, P(NN), Q(NN)
    COMPLEX(KIND=DWP), INTENT(INOUT) :: Z(LDZ,N)

    COMPLEX(KIND=DWP) :: T
    INTEGER :: IP, IQ, I

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(IP,IQ,I,T) SHARED(NN,Z,J,P,Q)
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

  SUBROUTINE ZSTEP_LOOP(N, U, LDU, A, LDA, Z, LDZ, J, SIGMA, NN, P, Q, R, NM, DZ, N_2, STEP, TT, W, INFO)
#ifdef ANIMATE
    USE, INTRINSIC :: ISO_C_BINDING
#endif
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, LDU, LDA, LDZ, J(N), NN, P(NN), Q(NN), NM, N_2, TT
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(LDA,N)
    COMPLEX(KIND=DWP), INTENT(OUT) :: U(LDU,N), Z(LDZ,N), W(2,3,N_2)
    REAL(KIND=DWP), INTENT(OUT) :: SIGMA(N)
    TYPE(ZPROC), INTENT(IN) :: R
    TYPE(AW), INTENT(OUT) :: DZ(NM)
    INTEGER, INTENT(OUT) :: STEP(N_2), INFO

    INTEGER :: S, SL
#ifdef ANIMATE
    CHARACTER(LEN=8,KIND=c_char), PARAMETER :: GNAME = c_char_'zjkGstp'//C_NULL_CHAR, FNAME = c_char_'zjkFstp'//C_NULL_CHAR
    INTEGER(KIND=c_int), PARAMETER :: ACT = 2, SX = ANIMATE, SY = ANIMATE, BPP = 8
    INTEGER(KIND=c_intptr_t) :: CTX
    INTEGER(KIND=c_int) :: NF
    INTEGER(KIND=c_size_t) :: LDF
    INTEGER(KIND=c_intptr_t), EXTERNAL :: PVN_CVIS_START
    INTEGER(KIND=c_int), EXTERNAL :: PVN_CVIS_FRAME, PVN_CVIS_STOP
    LDF = INT(N,c_size_t)
    NF = INT(N,c_int)
#endif

#ifdef USE_FAST
    INFO = 0
#else
    IF (N .LT. 3) THEN
       INFO = -1
    ELSE IF (LDU .LT. N) THEN
       INFO = -3
    ELSE IF (LDA .LT. N) THEN
       INFO = -5
    ELSE IF (LDZ .LT. N) THEN
       INFO = -7
    ELSE IF (NN .LT. 0) THEN
       INFO = -10
    ELSE IF (NM .LT. NN) THEN
       INFO = -14
    ELSE IF (N_2 .LT. 0) THEN
       INFO = -16
    ELSE IF (TT .LT. 0) THEN
       INFO = -18
    ELSE IF (TT .GT. NN) THEN
       INFO = -18
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN
#endif

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(S,INFO) SHARED(N,U)
    DO S = 1, N
       DO INFO = 1, S-1
          U(INFO,S) = Z_ZERO
       END DO
       U(S,S) = Z_ONE
       DO INFO = S+1, N
          U(INFO,S) = Z_ZERO
       END DO
    END DO
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(S,INFO) SHARED(N,Z)
    DO S = 1, N
       DO INFO = 1, S-1
          Z(INFO,S) = Z_ZERO
       END DO
       Z(S,S) = Z_ONE
       DO INFO = S+1, N
          Z(INFO,S) = Z_ZERO
       END DO
    END DO
    !$OMP END PARALLEL DO

    INFO = 0
#ifdef ANIMATE
    CTX = PVN_CVIS_START(NF, NF, ACT, GNAME, FNAME)
    IF (CTX .EQ. 0_c_intptr_t) THEN
       INFO = -20
       WRITE (ERROR_UNIT,'(A)') 'PVN_CVIS_START'
       FLUSH(ERROR_UNIT)
       RETURN
    END IF
#endif

#ifndef USE_FAST
    WRITE (ERROR_UNIT,'(A)') '"STEP","OLDLEN","BUILDs","TRANSFs","NEWLEN","NORMF2D","TRIG","HYP"'
    FLUSH(ERROR_UNIT)
#endif

    S = 0
    DO WHILE (S .GE. 0)
#ifdef ANIMATE
       NF = PVN_CVIS_FRAME(CTX, A, LDF)
       IF (NF .NE. 0_c_int) THEN
          WRITE (ERROR_UNIT,'(A,I11)') 'PVN_CVIS_FRAME:', NF
          FLUSH(ERROR_UNIT)
          RETURN
       END IF
#endif
       CALL ZSTEP_EXEC(S, N, U, LDU, A, LDA, Z, LDZ, J, SIGMA, NN, P, Q, R, NM, DZ, TT, N_2, STEP, SL, W, INFO)
       IF (INFO .LT. 0) EXIT
       IF (SL .LE. 0) EXIT
       S = S + 1
    END DO
    IF (SL .LT. 0) THEN
       INFO = SL
    ELSE IF (INFO .GE. 0) THEN
       INFO = S
    END IF

#ifdef ANIMATE
    NF = PVN_CVIS_STOP(CTX, SX, SY, BPP, GNAME, BPP, FNAME)
    IF (NF .NE. 0_c_int) THEN
       WRITE (ERROR_UNIT,'(A,I11)') 'PVN_CVIS_STOP:', NF
       FLUSH(ERROR_UNIT)
    END IF
    CTX = 0_c_intptr_t
#endif

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(S) SHARED(N,SIGMA,A)
    DO S = 1, N
       SIGMA(S) = REAL(A(S,S))
#ifndef USE_FAST
       ! for a simple calculation of ||off(A)||
       A(S,S) = Z_ZERO
#endif
    END DO
    !$OMP END PARALLEL DO
    CALL JZHJ(N, Z, LDZ, J, NN, P, Q)
  END SUBROUTINE ZSTEP_LOOP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE zstep
