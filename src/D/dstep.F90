MODULE dstep
  !$ USE omp_lib
  USE dtransf
  USE dtypes
  USE jstep
  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE REAL(KIND=DWP) FUNCTION DMAGF2T(N, P, Q, A, LDA, J)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, P, Q, LDA, J(N)
    REAL(KIND=DWP), INTENT(IN) :: A(LDA,N)

    SELECT CASE (J(P) + J(Q))
    CASE (-2)
       IF ((A(Q,P) .NE. D_ZERO) .OR. (A(P,Q) .NE. D_ZERO) .OR. &
            (SIGN(D_ONE, A(P,P)) .EQ. D_MONE) .OR. (SIGN(D_ONE, A(Q,Q)) .EQ. D_MONE) .OR. &
            (A(Q,Q) .LT. A(P,P))) THEN
          DMAGF2T = DASUM2(A(Q,P), A(P,Q))
       ELSE ! no transform
          DMAGF2T = QUIET_NAN((P - 1) * N + (Q - 1))
       END IF
    CASE (2)
       IF ((A(Q,P) .NE. D_ZERO) .OR. (A(P,Q) .NE. D_ZERO) .OR. &
            (SIGN(D_ONE, A(P,P)) .EQ. D_MONE) .OR. (SIGN(D_ONE, A(Q,Q)) .EQ. D_MONE) .OR. &
            (A(P,P) .LT. A(Q,Q))) THEN
          DMAGF2T = DASUM2(A(Q,P), A(P,Q))
       ELSE ! no transform
          DMAGF2T = QUIET_NAN((P - 1) * N + (Q - 1))
       END IF
    CASE DEFAULT ! invalid J
       DMAGF2T = QUIET_NAN((P - 1) * N + (Q - 1))
    END SELECT
  END FUNCTION DMAGF2T

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE REAL(KIND=DWP) FUNCTION DMAGF2H(N, P, Q, A, LDA, J)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, P, Q, LDA, J(N)
    REAL(KIND=DWP), INTENT(IN) :: A(LDA,N)

    REAL(KIND=DWP) :: A2(2,2), U2(2,2), Z2(2,2)
    INTEGER :: J2(2), INFO

    IF ((A(Q,P) .NE. D_ZERO) .OR. (A(P,Q) .NE. D_ZERO) .OR. &
         (SIGN(D_ONE, A(P,P)) .EQ. D_MONE) .OR. (SIGN(D_ONE, A(Q,Q)) .EQ. D_MONE)) THEN
       A2(1,1) = A(P,P)
       A2(2,1) = A(Q,P)
       A2(1,2) = A(P,Q)
       A2(2,2) = A(Q,Q)
       J2(1) = J(P)
       J2(2) = J(Q)
       CALL DHSVD2(A2, J2, U2, Z2, TH1FIX, INFO)
       IF (INFO .LT. 0) THEN
          DMAGF2H = QUIET_NAN((P - 1) * N + (Q - 1))
       ELSE IF (IAND(INFO, 1) .EQ. 0) THEN
          DMAGF2H = D_ZERO
       ELSE IF (IAND(INFO, 4) .EQ. 0) THEN
          DMAGF2H = DASUM2(A(Q,P), A(P,Q))
       ELSE IF (IAND(INFO, 8) .EQ. 0) THEN
          DMAGF2H = ABODNDF2(Z2, N, A(1,P), A(1,Q), D_ZERO, D_ZERO, P, Q)
       ELSE ! TH1FIX
          DMAGF2H = ABODNDF2(Z2, N, A(1,P), A(1,Q), A2(2,1), A2(1,2), P, Q)
       END IF
    ELSE ! no transform
       DMAGF2H = QUIET_NAN((P - 1) * N + (Q - 1))
    END IF
  END FUNCTION DMAGF2H

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DPROC_INIT(NT, ID_NCP, ID_TRU, R, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT
    INTEGER, INTENT(INOUT) :: ID_NCP, ID_TRU
    TYPE(DPROC), INTENT(OUT) :: R
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

#ifdef USE_SUN
    R%NCP = ID_NCP
#else
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
#endif

    IF (ID_TRU .EQ. 0) ID_TRU = 1
#ifdef USE_SUN
    R%TRU = ID_TRU
#else
    SELECT CASE (ID_TRU)
    CASE (1)
       R%TRU => TRU1
    CASE (2)
       R%TRU => TRU2
    CASE DEFAULT
       R%TRU => NULL()
       INFO = -3
    END SELECT
#endif

    IF (NT .GT. 1) THEN
#ifdef USE_SUN
       R%SRT = 1
#else
       R%SRT => AW_SRT1
#endif
    ELSE ! <= 1 (single-threaded)
#ifdef USE_SUN
       R%SRT = 2
#else
       R%SRT => AW_SRT2
#endif
    END IF
  END SUBROUTINE DPROC_INIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DSTEP_BUILD(S, N, A, LDA, J, NN, P, Q, R, NM, DZ, TT, N_2, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: S, N, LDA, J(N), NN, P(NN), Q(NN), NM, TT, N_2
    REAL(KIND=DWP), INTENT(IN) :: A(LDA,N)
    TYPE(DPROC), INTENT(IN) :: R
    TYPE(AW), INTENT(OUT) :: DZ(NM)
    INTEGER, INTENT(OUT) :: SL, STEP(N_2)
    INTEGER(KIND=DWP), INTENT(OUT) :: INFO

    TYPE(AW) :: X
    INTEGER(KIND=DWP) :: ID
    INTEGER :: IP, IQ, I
    INTEGER :: II, IT
    !$ SAVE :: II, IT

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
       DZ(I)%W = DMAGF2T(N, IP, IQ, A, LDA, J)
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
       DZ(I)%W = DMAGF2H(N, IP, IQ, A, LDA, J)
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
    I = OPEN_LOG('DSTEP_BUILD', S)
#endif
#ifdef USE_SUN
    SELECT CASE (R%SRT)
    CASE (1)
       CALL AW_SRT1(IT, NM, DZ, ID)
    CASE (2)
       CALL AW_SRT2(IT, NM, DZ, ID)
    CASE DEFAULT
       STOP 'R%SRT'
    END SELECT
#else
    CALL R%SRT(IT, NM, DZ, ID)
#endif
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

#ifdef USE_SUN
    SELECT CASE (R%NCP)
    CASE (1)
       CALL AW_NCP1(N, IT, NM, DZ, MIN(IT, N_2), SL, STEP, ID)
    CASE (2)
       CALL AW_NCP2(N, IT, NM, DZ, MIN(IT, N_2), SL, STEP, ID)
    CASE (3)
       CALL AW_NCP3(N, IT, NM, DZ, MIN(IT, N_2), SL, STEP, ID)
    CASE DEFAULT
       STOP 'R%NCP'
    END SELECT
#else
    CALL R%NCP(N, IT, NM, DZ, MIN(IT, N_2), SL, STEP, ID)
#endif
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
  END SUBROUTINE DSTEP_BUILD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DSTEP_TRANSF(N, U, LDU, A, LDA, Z, LDZ, J, SIGMA, NN, DZ, SL, STEP, W, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, LDU, LDA, LDZ, J(N), NN, SL, STEP(SL)
    REAL(KIND=DWP), INTENT(INOUT) :: U(LDU,N), A(LDA,N), Z(LDZ,N)
    REAL(KIND=DWP), INTENT(OUT) :: SIGMA(N), W(2,3,SL)
    TYPE(AW), INTENT(INOUT) :: DZ(NN)
    INTEGER, INTENT(OUT) :: INFO

    REAL(KIND=DWP) :: V(2,2), B(2,2)
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

       CALL DHSVD2(B, K, V, W(:,:,I), TH1FIX, M)
       DZ(STEP(I))%I = M
       IF (M .GE. 0) THEN
          W(1,3,I) = B(1,1)
          W(2,3,I) = B(2,2)
          IF (IAND(M, 2) .NE. 0) THEN
             CALL BA(V, N, A(P,1), A(Q,1), LDA)
             CALL UT2(V)
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
             A(Q,P) = D_ZERO
             A(P,Q) = D_ZERO
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
  END SUBROUTINE DSTEP_TRANSF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DSTEP_EXEC(S, N, U, LDU, A, LDA, Z, LDZ, J, SIGMA, NN, P, Q, R, NM, DZ, TT, N_2, STEP, SL, W, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: S, N, LDU, LDA, LDZ, J(N), NN, P(NN), Q(NN), NM, TT, N_2
    REAL(KIND=DWP), INTENT(INOUT) :: U(LDU,N), A(LDA,N), Z(LDZ,N)
    REAL(KIND=DWP), INTENT(OUT) :: SIGMA(N), W(2,3,N_2)
    TYPE(DPROC), INTENT(IN) :: R
    TYPE(AW), INTENT(OUT) :: DZ(NM)
    INTEGER, INTENT(OUT) :: STEP(N_2), SL, INFO

    INTEGER(KIND=DWP) :: IB, IT
    INTEGER :: NL

#ifndef USE_FAST
    WRITE (ERROR_UNIT,'(I11,A)',ADVANCE='NO') S, ','
    FLUSH(ERROR_UNIT)
#endif

    SL = 0
    CALL DSTEP_BUILD(S, N, A, LDA, J, NN, P, Q, R, NM, DZ, TT, N_2, SL, STEP, IB)
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
    CALL DSTEP_TRANSF(N, U, LDU, A, LDA, Z, LDZ, J, SIGMA, NN, DZ, SL, STEP, W, NL)
#ifndef USE_FAST
    IT = GET_SYS_TLAP(IT)
    WRITE (ERROR_UNIT,'(F12.6,A,I11,A,ES25.17E3,2(A,I11))') &
         (IT / REAL(GET_SYS_TRES(),DWP)), ',', NL, ',', SIGMA(1), ',',INT(SIGMA(2)),',',INT(SIGMA(3))
    FLUSH(ERROR_UNIT)
#endif

    SL = MIN(SL, NL)
  END SUBROUTINE DSTEP_EXEC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE JZTJ(N, Z, LDZ, J, NN, P, Q)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, LDZ, J(N), NN, P(NN), Q(NN)
    REAL(KIND=DWP), INTENT(INOUT) :: Z(LDZ,N)

    REAL(KIND=DWP) :: T
    INTEGER :: IP, IQ, I

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(IP,IQ,I,T) SHARED(NN,Z,J,P,Q)
    DO I = 1, NN
       IP = P(I)
       IQ = Q(I)
       T = Z(IP,IQ)
       Z(IP,IQ) = Z(IQ,IP)
       Z(IQ,IP) = T
       IF (J(IP) .NE. J(IQ)) THEN
          Z(IP,IQ) = -Z(IP,IQ)
          Z(IQ,IP) = -Z(IQ,IP)
       END IF
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE JZTJ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DSTEP_LOOP(N, U, LDU, A, LDA, Z, LDZ, J, SIGMA, NN, P, Q, R, NM, DZ, N_2, STEP, TT, W, INFO)
#ifdef ANIMATE
    USE vn_mtxvis_f
#endif
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, LDU, LDA, LDZ, J(N), NN, P(NN), Q(NN), NM, N_2, TT
    REAL(KIND=DWP), INTENT(INOUT) :: A(LDA,N)
    REAL(KIND=DWP), INTENT(OUT) :: U(LDU,N), Z(LDZ,N), SIGMA(N), W(2,3,N_2)
    TYPE(DPROC), INTENT(IN) :: R
    TYPE(AW), INTENT(OUT) :: DZ(NM)
    INTEGER, INTENT(OUT) :: STEP(N_2), INFO

    INTEGER :: S, SL
#ifdef ANIMATE
    CHARACTER(LEN=7), PARAMETER :: FNAME = 'djkGstp'
    INTEGER, PARAMETER :: ACT = IOR(IOR(VN_MTXVIS_OP_A, VN_MTXVIS_FN_Lg), VN_MTXVIS_FF_Bin)
    INTEGER, PARAMETER :: SX = ANIMATE, SY = ANIMATE
    TYPE(c_ptr) :: CTX
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
          U(INFO,S) = D_ZERO
       END DO
       U(S,S) = D_ONE
       DO INFO = S+1, N
          U(INFO,S) = D_ZERO
       END DO
    END DO
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(S,INFO) SHARED(N,Z)
    DO S = 1, N
       DO INFO = 1, S-1
          Z(INFO,S) = D_ZERO
       END DO
       Z(S,S) = D_ONE
       DO INFO = S+1, N
          Z(INFO,S) = D_ZERO
       END DO
    END DO
    !$OMP END PARALLEL DO

#ifdef ANIMATE
    INFO = VN_MTXVIS_START(CTX, FNAME, ACT, N, N, SX, SY, LEN_TRIM(FNAME))
    IF (INFO .NE. 0) THEN
       WRITE (ERROR_UNIT,'(A,I11)') 'VN_MTXVIS_START:', INFO
       FLUSH(ERROR_UNIT)
       RETURN
    END IF
#else
    INFO = 0
#endif

#ifndef USE_FAST
    WRITE (ERROR_UNIT,'(A)') '"STEP","OLDLEN","BUILDs","TRANSFs","NEWLEN","NORMF2D","TRIG","HYP"'
    FLUSH(ERROR_UNIT)
#endif

    S = 0
    DO WHILE (S .GE. 0)
#ifdef ANIMATE
       INFO = VN_MTXVIS_FRAME(CTX, A, N)
       IF (INFO .NE. 0) THEN
          WRITE (ERROR_UNIT,'(A,I11)') 'VN_MTXVIS_FRAME:', INFO
          FLUSH(ERROR_UNIT)
          RETURN
       END IF
#endif
       CALL DSTEP_EXEC(S, N, U, LDU, A, LDA, Z, LDZ, J, SIGMA, NN, P, Q, R, NM, DZ, TT, N_2, STEP, SL, W, INFO)
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
    S = VN_MTXVIS_STOP(CTX)
    IF (S .NE. 0) THEN
       WRITE (ERROR_UNIT,'(A,I11)') 'VN_MTXVIS_STOP:', S
       FLUSH(ERROR_UNIT)
    END IF
#endif

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(S) SHARED(N,SIGMA,A)
    DO S = 1, N
       SIGMA(S) = A(S,S)
#ifndef USE_FAST
       ! for a simple calculation of ||off(A)||
       A(S,S) = D_ZERO
#endif
    END DO
    !$OMP END PARALLEL DO
    CALL JZTJ(N, Z, LDZ, J, NN, P, Q)
  END SUBROUTINE DSTEP_LOOP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE dstep
