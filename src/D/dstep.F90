MODULE DSTEP
  USE, INTRINSIC :: ISO_C_BINDING
  USE OMP_LIB
  USE DTRANSF
  USE DTYPES
  USE JSTEP
  USE UTILS
  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE REAL(KIND=DWP) FUNCTION DMAGF2T(N, P, Q, A, LDA, J)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, P, Q, LDA, J(N)
    REAL(KIND=DWP), INTENT(IN) :: A(LDA,N)

    REAL(KIND=DWP) :: A2(2,2), U2(2,2), Z2(2,2)
    INTEGER :: INFO

    INFO = J(P) + J(Q)
    IF (INFO .EQ. 2) THEN
       IF ((A(Q,P) .NE. D_ZERO) .OR. (A(P,Q) .NE. D_ZERO) .OR. (SIGN(D_ONE, A(P,P)) .EQ. D_MONE) .OR. &
            (SIGN(D_ONE, A(Q,Q)) .EQ. D_MONE) .OR. (A(P,P) .LT. A(Q,Q))) THEN
          A2(2,1) = ABS(A(Q,P))
          A2(1,2) = ABS(A(P,Q))
          IF (A2(2,1) .GE. A2(1,2)) THEN
             IF (A2(2,1) .EQ. D_ZERO) THEN
                DMAGF2T = D_ZERO
             ELSE ! A2(2,1) .GT. D_ZERO
                A2(1,1) = A2(1,2) / A2(2,1)
                DMAGF2T = A2(1,1) * A2(1,2) + A2(2,1)
                DMAGF2T = DMAGF2T * A2(2,1)
             END IF
          ELSE ! A2(2,1) .LT. A2(1,2)
             A2(2,2) = A2(2,1) / A2(1,2)
             DMAGF2T = A2(1,2) + A2(2,1) * A2(2,2)
             DMAGF2T = DMAGF2T * A2(1,2)
          END IF
       ELSE ! no transform
          DMAGF2T = QUIET_NAN((P - 1) * N + (Q - 1))
       END IF
    ELSE IF (INFO .EQ. -2) THEN
       IF ((A(Q,P) .NE. D_ZERO) .OR. (A(P,Q) .NE. D_ZERO) .OR. (SIGN(D_ONE, A(P,P)) .EQ. D_MONE) .OR. &
            (SIGN(D_ONE, A(Q,Q)) .EQ. D_MONE) .OR. (A(Q,Q) .LT. A(P,P))) THEN
          A2(2,1) = ABS(A(Q,P))
          A2(1,2) = ABS(A(P,Q))
          IF (A2(2,1) .GE. A2(1,2)) THEN
             IF (A2(2,1) .EQ. D_ZERO) THEN
                DMAGF2T = D_ZERO
             ELSE ! A2(2,1) .GT. D_ZERO
                A2(1,1) = A2(1,2) / A2(2,1)
                DMAGF2T = A2(1,1) * A2(1,2) + A2(2,1)
                DMAGF2T = DMAGF2T * A2(2,1)
             END IF
          ELSE ! A2(2,1) .LT. A2(1,2)
             A2(2,2) = A2(2,1) / A2(1,2)
             DMAGF2T = A2(1,2) + A2(2,1) * A2(2,2)
             DMAGF2T = DMAGF2T * A2(1,2)
          END IF
       ELSE ! no transform
          DMAGF2T = QUIET_NAN((P - 1) * N + (Q - 1))
       END IF
    ELSE ! invalid J
       DMAGF2T = QUIET_NAN((P - 1) * N + (Q - 1))
    END IF
  END FUNCTION DMAGF2T

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! #ifdef USE_EXTENDED
!   REAL(KIND=DWP) FUNCTION DMAGF2H(N, P, Q, A, LDA, J)
! #else
  PURE REAL(KIND=DWP) FUNCTION DMAGF2H(N, P, Q, A, LDA, J)
! #endif
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, P, Q, LDA, J(N)
    REAL(KIND=DWP), INTENT(IN) :: A(LDA,N)

    REAL(KIND=DWP) :: A2(2,2), U2(2,2), Z2(2,2)
    INTEGER :: J2(2), INFO

! #ifdef USE_EXTENDED
!     EXTERNAL :: DHSVD2
! #endif

    IF ((A(Q,P) .NE. D_ZERO) .OR. (A(P,Q) .NE. D_ZERO) .OR. (SIGN(D_ONE, A(P,P)) .EQ. D_MONE) .OR. &
         (SIGN(D_ONE, A(Q,Q)) .EQ. D_MONE)) THEN
       A2(1,1) = A(P,P)
       A2(2,1) = A(Q,P)
       A2(1,2) = A(P,Q)
       A2(2,2) = A(Q,Q)
       J2(1) = J(P)
       J2(2) = J(Q)
       CALL DHSVD2(A2, J2, U2, Z2, INFO)
       IF (INFO .LE. 1) THEN
          DMAGF2H = QUIET_NAN((P - 1) * N + (Q - 1))
       ELSE IF (IAND(INFO, 1) .EQ. 0) THEN
          A2(2,1) = ABS(A(Q,P))
          A2(1,2) = ABS(A(P,Q))
          IF (A2(2,1) .GE. A2(1,2)) THEN
             IF (A2(2,1) .EQ. D_ZERO) THEN
                DMAGF2H = D_ZERO
             ELSE ! A2(2,1) .GT. D_ZERO
                A2(1,1) = A2(1,2) / A2(2,1)
                DMAGF2H = A2(1,1) * A2(1,2) + A2(2,1)
                DMAGF2H = DMAGF2H * A2(2,1)
             END IF
          ELSE ! A2(2,1) .LT. A2(1,2)
             A2(2,2) = A2(2,1) / A2(1,2)
             DMAGF2H = A2(1,2) + A2(2,1) * A2(2,2)
             DMAGF2H = DMAGF2H * A2(1,2)
          END IF
       ELSE ! a non-trivial transform
          DMAGF2H = ABODNDF2(Z2, N, A(1,P), A(1,Q), P, Q)
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
  END SUBROUTINE DPROC_INIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DSTEP_BUILD(NT, S, N, A, LDA, J, NN, P, Q, R, NM, DZ, TT, N_2, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, S, N, LDA, J(N), NN, P(NN), Q(NN), NM, TT, N_2
    REAL(KIND=DWP), INTENT(IN) :: A(LDA,N)
    TYPE(DPROC), INTENT(IN) :: R
    TYPE(AW), INTENT(OUT), TARGET :: DZ(NM)
    INTEGER, INTENT(OUT) :: SL, STEP(N_2), INFO

    TYPE(AW) :: X
    INTEGER :: IP, IQ, I, II, IT
#ifndef NDEBUG
    REAL(KIND=DWP) :: T
#endif

    SL = 0
#ifndef NDEBUG
    I = -1
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
    ELSE IF (TT .LT. 0) THEN
       INFO = -13
    ELSE IF (TT .GT. NN) THEN
       INFO = -13
    ELSE IF (N_2 .LT. 0) THEN
       INFO = -14
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    INFO = GET_SYS_US()
    IF (N .EQ. 0) GOTO 1
    IF (NN .EQ. 0) GOTO 1
    IF (N_2 .EQ. 0) GOTO 1

    IT = 0
    !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(IP,IQ,I) SHARED(TT,N,A,LDA,J,P,Q,R,DZ) REDUCTION(+:IT)
    DO I = 1, TT
       IP = P(I)
       IQ = Q(I)
       DZ(I)%P = IP
       DZ(I)%Q = IQ
       DZ(I)%B = IQ - IP
       DZ(I)%W = DMAGF2T(N, IP, IQ, A, LDA, J)
       IF (DZ(I)%W .EQ. DZ(I)%W) IT = IT + 1
    END DO
    !$OMP END PARALLEL DO

    II = 0
    !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(IP,IQ,I) SHARED(NN,TT,N,A,LDA,J,P,Q,R,DZ) REDUCTION(+:II)
    DO I = TT+1, NN
       IP = P(I)
       IQ = Q(I)
       DZ(I)%P = IP
       DZ(I)%Q = IQ
       DZ(I)%B = IQ - IP
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
       DO WHILE (.NOT. (DZ(II)%W .EQ. DZ(II)%W))
          II = II - 1
       END DO
       DO WHILE (I .LT. II)
          IF (DZ(I)%W .EQ. DZ(I)%W) THEN
             I = I + 1
          ELSE ! NaN
             X = DZ(I)
             DZ(I) = DZ(II)
             DZ(II) = X
             I = I + 1
             II = II - 1
             DO WHILE (.NOT. (DZ(II)%W .EQ. DZ(II)%W))
                II = II - 1
             END DO
          END IF
       END DO
    END IF
#ifndef NDEBUG
    I = OPEN_LOG('DSTEP_BUILD', S)
#endif
    CALL R%SRT(NT, IT, NM, DZ, AW_CMP, II)
    IF (II .LT. 0) THEN
       INFO = -11
       RETURN
    END IF
#ifndef NDEBUG
    IF (I .NE. -1) THEN
       T = II * DUS2S
       CALL AW_OUT(I, '', IT, DZ, 0, STEP, II)
       WRITE (I,'(A,F12.6,A)',ADVANCE='NO') 'SORT: ', T, ' s, '
    END IF
#endif

    CALL R%NCP(NT, N, IT, NM, DZ, MIN(IT, N_2), SL, STEP, II)
    IF (II .LT. 0) THEN
       INFO = -13
       RETURN
    END IF
#ifndef NDEBUG
    IF (I .NE. -1) THEN
       T = II * DUS2S
       WRITE (I,'(A,F12.6,A)',ADVANCE='NO') 'NCP: ', T, ' s, '
    END IF
#endif

1   INFO = MAX(GET_SYS_US() - INFO, 1)
#ifndef NDEBUG
    IF (I .NE. -1) THEN
       T = INFO * DUS2S
       WRITE (I,'(A,F12.6,A)') 'BUILD: ', T, ' s'
       CALL AW_OUT(I, '', IT, DZ, SL, STEP, II)
       CLOSE(UNIT=I, IOSTAT=II)
    END IF
#endif
  END SUBROUTINE DSTEP_BUILD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DSTEP_TRANSF(NT, N, U, LDU, A, LDA, Z, LDZ, J, SIGMA, NN, NM, DZ, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, N, LDU, LDA, LDZ, J(N), NN, NM, SL, STEP(SL)
    REAL(KIND=DWP), INTENT(INOUT) :: U(LDU,N), A(LDA,N), Z(LDZ,N)
    REAL(KIND=DWP), INTENT(OUT), TARGET :: SIGMA(N)
    TYPE(AW), INTENT(INOUT), TARGET :: DZ(NM)
    INTEGER, INTENT(OUT) :: INFO

    REAL(KIND=DWP) :: V(2,2), B(2,2)
    REAL(KIND=DWP), POINTER, CONTIGUOUS :: W(:,:,:)
    REAL(KIND=DWP) :: TW
    INTEGER :: I, P, Q, K(2)
    INTEGER, POINTER, CONTIGUOUS :: IT(:)

! #ifdef USE_EXTENDED
!     EXTERNAL :: DHSVD2
! #endif

    CALL C_F_POINTER(C_LOC(DZ(1+NN)), W, [2,3,SL])
    CALL C_F_POINTER(C_LOC(SIGMA(1+N/2)), IT, [SL])
    INFO = HUGE(INFO)

    !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(I,P,Q,V,B,K) SHARED(N,SL,STEP,IT,DZ,J,U,A,Z,LDA,W) &
    !$OMP& REDUCTION(MIN:INFO)
    DO I = 1, SL
       P = DZ(STEP(I))%P
       Q = DZ(STEP(I))%Q

       B(1,1) = A(P,P)
       B(2,1) = A(Q,P)
       B(1,2) = A(P,Q)
       B(2,2) = A(Q,Q)
       K(1) = J(P)
       K(2) = J(Q)

       CALL DHSVD2(B, K, V, W(:,:,I), IT(I))
       INFO = IT(I)
       IF (IT(I) .GE. 0) THEN
          W(1,3,I) = B(1,1)
          W(2,3,I) = B(2,2)
          IF (IAND(IT(I), 2) .NE. 0) THEN
             CALL BA(V, N, A(P,1), A(Q,1), LDA)
             CALL UT2(V)
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
          IF ((IAND(IT(I), 1) .NE. 0) .AND. ((IAND(IT(I), 2) .NE. 0) .OR. (IAND(IT(I), 4) .NE. 0))) INFO = INFO + 1

          A(P,P) = W(1,3,I)
          A(Q,P) = D_ZERO
          A(P,Q) = D_ZERO
          A(Q,Q) = W(2,3,I)
       END DO
       !$OMP END PARALLEL DO
    END IF

    TW = D_MZERO
    P = 0
    Q = 0
    IF (INFO .GE. 0) THEN
       !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(I) SHARED(SL,STEP,DZ,IT,J) REDUCTION(+:TW,P,Q)
       DO I = 1, SL
          IF ((IAND(IT(I), 2) .NE. 0) .OR. (IAND(IT(I), 4) .NE. 0)) THEN
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
  END SUBROUTINE DSTEP_TRANSF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DSTEP_EXEC(NT, S, N, U, LDU, A, LDA, Z, LDZ, J, SIGMA, NN, P, Q, R, NM, DZ, TT, N_2, STEP, SL, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, S, N, LDU, LDA, LDZ, J(N), NN, P(NN), Q(NN), NM, TT, N_2
    REAL(KIND=DWP), INTENT(INOUT) :: U(LDU,N), A(LDA,N), Z(LDZ,N)
    REAL(KIND=DWP), INTENT(OUT), TARGET :: SIGMA(N)
    TYPE(DPROC), INTENT(IN) :: R
    TYPE(AW), INTENT(OUT), TARGET :: DZ(NM)
    INTEGER, INTENT(OUT) :: STEP(N_2), SL, INFO

    INTEGER :: IT, NL

    SL = 0
    INFO = 0
#ifndef NDEBUG
    IF (IsCtrlC()) RETURN
#endif
    WRITE (ERROR_UNIT,'(I10,A)',ADVANCE='NO') S, ','
    FLUSH(ERROR_UNIT)
    CALL DSTEP_BUILD(NT, S, N, A, LDA, J, NN, P, Q, R, NM, DZ, TT, N_2, SL, STEP, INFO)
    IF (INFO .LE. 0) THEN
       IT = INFO
    ELSE
       IT = SL
    END IF
    WRITE (ERROR_UNIT,'(I11,A)',ADVANCE='NO') IT, ','
    FLUSH(ERROR_UNIT)

    IF (INFO .LE. 0) THEN
       WRITE (ERROR_UNIT,'(F12.6,A,F12.6,A,I11,A,ES25.17E3,2(A,I11))') &
            D_MZERO, ',', D_ZERO, ',', -1, ',', D_MZERO, ',',0,',',0
    ELSE IF (SL .LE. 0) THEN
       WRITE (ERROR_UNIT,'(F12.6,A,F12.6,A,I11,A,ES25.17E3,2(A,I11))') &
            (INFO * DUS2S), ',', D_ZERO, ',', -1, ',', D_MZERO, ',',0,',',0
    ELSE ! default
       WRITE (ERROR_UNIT,'(F12.6,A)',ADVANCE='NO') (INFO * DUS2S), ','
    END IF
    FLUSH(ERROR_UNIT)

    IF (SL .LE. 0) RETURN
#ifndef NDEBUG
    IF (IsCtrlC()) INFO = 0
#endif
    IF (INFO .LE. 0) THEN
       WRITE (ERROR_UNIT,'(F12.6,A,I11,A,ES25.17E3,2(A,I11))') &
            D_MZERO, ',', -1, ',', D_MZERO, ',',0,',',0
       FLUSH(ERROR_UNIT)
       RETURN
    END IF

    IT = GET_SYS_US()
    CALL DSTEP_TRANSF(NT, N, U, LDU, A, LDA, Z, LDZ, J, SIGMA, NN, NM, DZ, SL, STEP, NL)
    IT = MAX(GET_SYS_US() - IT, 1)
    WRITE (ERROR_UNIT,'(F12.6,A,I11,A,ES25.17E3,2(A,I11))') &
         (IT * DUS2S), ',', NL, ',', SIGMA(1), ',',INT(SIGMA(2)),',',INT(SIGMA(3))
    FLUSH(ERROR_UNIT)

    SL = MIN(SL, NL)
    INFO = INFO + IT
  END SUBROUTINE DSTEP_EXEC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE JZTJ(NT, N, Z, LDZ, J, NN, P, Q)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, N, LDZ, J(N), NN, P(NN), Q(NN)
    REAL(KIND=DWP), INTENT(INOUT) :: Z(LDZ,N)

    REAL(KIND=DWP) :: T
    INTEGER :: IP, IQ, I

    !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(IP,IQ,I,T) SHARED(NN,Z,J,P,Q)
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

  SUBROUTINE DSTEP_LOOP(NT, N, U, LDU, A, LDA, Z, LDZ, J, SIGMA, NN, P, Q, R, NM, DZ, N_2, STEP, TT, INFO)
#ifdef ANIMATE
    USE VN_MTXVIS_F
#endif
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, N, LDU, LDA, LDZ, J(N), NN, P(NN), Q(NN), NM, N_2, TT
    REAL(KIND=DWP), INTENT(INOUT) :: A(LDA,N)
    REAL(KIND=DWP), INTENT(OUT), TARGET :: U(LDU,N), Z(LDZ,N), SIGMA(N)
    TYPE(DPROC), INTENT(IN) :: R
    TYPE(AW), INTENT(OUT), TARGET :: DZ(NM)
    INTEGER, INTENT(OUT) :: STEP(N_2), INFO

    INTEGER :: S, SL
#ifdef ANIMATE
    CHARACTER(LEN=7), PARAMETER :: FNAME = 'djkGstp'
    INTEGER, PARAMETER :: ACT = IOR(IOR(VN_MTXVIS_OP_A, VN_MTXVIS_FN_Lg), VN_MTXVIS_FF_Bin)
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
    ELSE IF (TT .LT. 0) THEN
       INFO = -19
    ELSE IF (TT .GT. NN) THEN
       INFO = -19
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(S,INFO) SHARED(N,U)
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

    !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(S,INFO) SHARED(N,Z)
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

    WRITE (ERROR_UNIT,'(A)') '"STEP","OLDLEN","BUILDs","TRANSFs","NEWLEN","NORMF2D","TRIG","HYP"'
    FLUSH(ERROR_UNIT)

    S = 0
#ifdef NDEBUG
    DO WHILE (S .GE. 0)
#else
    DO WHILE ((S .GE. 0) .AND. (.NOT. IsCtrlC()))
#endif
#ifdef ANIMATE
       INFO = VN_MTXVIS_FRAME(CTX, A, N)
       IF (INFO .NE. 0) THEN
          WRITE (ERROR_UNIT,'(A,I11)') 'VN_MTXVIS_FRAME:', INFO
          FLUSH(ERROR_UNIT)
          RETURN
       END IF
#endif
       CALL DSTEP_EXEC(NT, S, N, U, LDU, A, LDA, Z, LDZ, J, SIGMA, NN, P, Q, R, NM, DZ, TT, N_2, STEP, SL, INFO)
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
    S = VN_MTXVIS_STOP(CTX)
    IF (S .NE. 0) THEN
       WRITE (ERROR_UNIT,'(A,I11)') 'VN_MTXVIS_STOP:', S
       FLUSH(ERROR_UNIT)
    END IF
#endif

    !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(S) SHARED(N,SIGMA,A)
    DO S = 1, N
       SIGMA(S) = A(S,S)
       ! for a simple calculation of ||off(A)||
       A(S,S) = D_ZERO
    END DO
    !$OMP END PARALLEL DO
    CALL JZTJ(NT, N, Z, LDZ, J, NN, P, Q)
  END SUBROUTINE DSTEP_LOOP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE DSTEP
