MODULE atypes
  !$ USE omp_lib
#ifndef USE_FAST
  USE timer
#endif
  USE utils

  IMPLICIT NONE

  ABSTRACT INTERFACE
     PURE SUBROUTINE ATRU(N, J, NN, P, Q, PN, INFO)
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: N, J(N), NN
       INTEGER, INTENT(OUT) :: P(NN), Q(NN), PN(2), INFO
     END SUBROUTINE ATRU
  END INTERFACE

  TYPE :: AW
     REAL(KIND=DWP) :: W ! weight
     INTEGER :: P ! row
     INTEGER :: Q ! column
     INTEGER :: B ! band (Q - P) > 0
     INTEGER :: I ! info
  END TYPE AW

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef NDEBUG
  SUBROUTINE AW_OUT(OU, HDR, NN, DZ, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: OU, NN, SL, STEP(SL)
    CHARACTER(LEN=*), INTENT(IN) :: HDR
    TYPE(AW), INTENT(IN) :: DZ(NN)
    INTEGER, INTENT(OUT) :: INFO

    INTEGER :: I, J

    IF (OU .EQ. -1) THEN
       INFO = -1
    ELSE IF (NN .LT. 0) THEN
       INFO = -3
    ELSE IF (SL .LT. 0) THEN
       INFO = -5
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

#ifndef USE_FAST
    INFO = GET_SYS_TIME()
#endif

    IF (LEN_TRIM(HDR) .GT. 0) THEN
       WRITE (OU,'(A)') TRIM(HDR)
    ELSE IF (SL .EQ. 0) THEN
       WRITE (OU,'(A)') '"J","W","P","Q","B"'
    ELSE ! SL > 0
       WRITE (OU,'(A)') '"I","J","W","P","Q","B"'
    END IF

    IF (SL .EQ. 0) THEN
       DO J = 1, NN
          WRITE (OU,'(I10,A,ES25.17E3,3(A,I10))') J, ',', DZ(J)%W, ',', DZ(J)%P, ',', DZ(J)%Q, ',', DZ(J)%B
       END DO
    ELSE ! SL > 0
       DO I = 1, SL
          J = STEP(I)
          WRITE (OU,'(2(I10,A),ES25.17E3,3(A,I10))') I, ',', J, ',', DZ(J)%W, ',', DZ(J)%P, ',', DZ(J)%Q, ',', DZ(J)%B
       END DO
    END IF

#ifndef USE_FAST
    INFO = GET_SYS_TIME() - INFO
#endif
  END SUBROUTINE AW_OUT
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! sorting by magnitude
  !   same magnitude ==> sorting by subdiagonals (bands)
  !     outer bands first
  !       within a band, the lower elements first
  PURE INTEGER FUNCTION AW_CMP(A, B)
    IMPLICIT NONE
    TYPE(AW), INTENT(IN) :: A, B

    AW_CMP = 0

    IF (A%W .LT. B%W) THEN
       AW_CMP = 1
    ELSE IF (A%W .GT. B%W) THEN
       AW_CMP = -1
    ELSE IF (A%W .EQ. B%W) THEN
       IF (A%B .LT. B%B) THEN
          AW_CMP = 2
       ELSE IF (A%B .GT. B%B) THEN
          AW_CMP = -2
       ELSE ! A%B = B%B
          IF (A%P .LT. B%P) THEN
             AW_CMP = 3
          ELSE IF (A%P .GT. B%P) THEN
             AW_CMP = -3
          ELSE ! A%P = B%P
             IF (A%Q .LT. B%Q) THEN
                AW_CMP = 4
             ELSE IF (A%Q .GT. B%Q) THEN
                AW_CMP = -4
             ELSE ! A%Q = B%Q
                AW_CMP = 0
             END IF
          END IF
       END IF
    ELSE ! NaN magnitude(s)
       IF (B%W .EQ. B%W) THEN
          AW_CMP = 5
       ELSE ! NaN(B%W)
          AW_CMP = -5
       END IF
    END IF
  END FUNCTION AW_CMP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE AW_MERGE(N, M, A, I, J, B, K, L)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, M, I, J, K
    TYPE(AW), INTENT(IN) :: A(*)
    TYPE(AW), INTENT(OUT) :: B(*)
    INTEGER, INTENT(OUT) :: L

    INTEGER :: II, JJ, KK, IL, JL, KL, LL

    L = 0
    IL = I + N - 1
    JL = J + M - 1
    KL = K + N + M - 1

    II = I
    JJ = J
    KK = K

    DO WHILE (KK .LE. KL)
       IF (II .GT. IL) THEN
          DO WHILE (KK .LE. KL)
             B(KK) = A(JJ)
             JJ = JJ + 1
             KK = KK + 1
          END DO
          EXIT
       END IF
       IF (JJ .GT. JL) THEN
          DO WHILE (KK .LE. KL)
             B(KK) = A(II)
             II = II + 1
             KK = KK + 1
          END DO
          EXIT
       END IF
       LL = AW_CMP(A(II), A(JJ))
       IF (LL .LT. 0) THEN
          B(KK) = A(II)
          II = II + 1
       ELSE IF (LL .EQ. 0) THEN
          B(KK) = A(II)
          II = II + 1
          KK = KK + 1
          B(KK) = A(JJ)
          JJ = JJ + 1
       ELSE ! .GT. 0
          B(KK) = A(JJ)
          JJ = JJ + 1
          L = L + 1
       END IF
       KK = KK + 1
    END DO
  END SUBROUTINE AW_MERGE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE AW_SORT(N, A, B)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    TYPE(AW), INTENT(INOUT) :: A(N)
    TYPE(AW), INTENT(OUT) :: B(N)

    INTEGER :: I, J, M, W, W2
    LOGICAL :: D

    W = 1
    D = .FALSE.
    DO WHILE (W .LT. N)
       M = W
       W2 = W * 2
       IF (D) THEN
          DO I = 1, N, W2
             J = I + M
             CALL AW_MERGE(MIN(M, N+1-I), MAX(MIN(M, N+1-J), 0), B, I, J, A, I, W)
          END DO
       ELSE ! A -> B
          DO I = 1, N, W2
             J = I + M
             CALL AW_MERGE(MIN(M, N+1-I), MAX(MIN(M, N+1-J), 0), A, I, J, B, I, W)
          END DO
       END IF
       W = W2
       D = (.NOT. D)
    END DO

    IF (D) THEN
       DO I = 1, N
          A(I) = B(I)
       END DO
    END IF
  END SUBROUTINE AW_SORT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE AW_SRT1(NN, NM, DZ, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NN, NM
    TYPE(AW), INTENT(INOUT) :: DZ(NM)
    INTEGER, INTENT(OUT) :: INFO

    INTEGER :: NT, T, I, J, K, TE
    ! elements per array, thread, # of + comparisons
    ! A => DZ(1:EPA), B => DZ(EPA+1:NM)
    INTEGER :: EPA, EPT, L
    !$ SAVE :: EPA, EPT, L

    IF (NN .LT. 0) THEN
       INFO = -1
    ELSE IF (NM .LT. NN) THEN
       INFO = -2
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN
    IF (NN .EQ. 0) RETURN

#ifndef USE_FAST
    INFO = GET_SYS_TIME()
#endif

    NT = 1
    !$ NT = MAX(NT, OMP_GET_MAX_THREADS())
    I = MOD(NN, NT)
    IF (I .GT. 0) I = NT - I
    EPA = NN + I
    EPT = EPA / NT

    ! virtual elements
    DO I = NN+1, EPA
       DZ(I)%W = QUIET_NAN(I)
       DZ(I)%P = 0
       DZ(I)%Q = 0
       DZ(I)%B = 0
       DZ(I)%I = 0
    END DO

    ! Baudet-Stevenson odd-even sort with merge-splitting of the subarrays; see:
    ! Baudet and Stevenson, Optimal Sorting Algorithms for Parallel Computers,
    ! IEEE Transactions on Computers, C-27(1):84--87, Jan 1978.
    ! doi:10.1109/TC.1978.1674957

    !$OMP PARALLEL NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(I) SHARED(EPA,EPT,DZ)
    I = 1
    !$ I = INT(OMP_GET_THREAD_NUM()) * EPT + I
    CALL AW_SORT(EPT, DZ(I), DZ(I+EPA))
    !$OMP END PARALLEL

    DO WHILE (.TRUE.)
       TE = 0

       L = 0
       !$OMP PARALLEL NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(T,I,J,K) SHARED(NT,EPA,EPT,DZ) REDUCTION(+:L)
       T = 0
       !$ T = INT(OMP_GET_THREAD_NUM())
       I = T * EPT + 1
       K = I + EPA
       IF (MOD(T, 2) .EQ. 0) THEN
          ! merge with T + 1
          J = I + EPT
          CALL AW_MERGE(EPT, MIN(EPT, EPA+1-J), DZ, I, J, DZ, K, L)
       ELSE ! nothing to do
          L = 0
       END IF
       !$OMP END PARALLEL
       TE = TE + L

       L = 0
       !$OMP PARALLEL NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(T,I,J,K) SHARED(NT,EPA,EPT,DZ) REDUCTION(+:L)
       T = 0
       !$ T = INT(OMP_GET_THREAD_NUM())
       K = T * EPT + 1
       I = K + EPA
       IF (MOD(T, 2) .EQ. 1) THEN
          ! merge with T + 1
          J = I + EPT
          CALL AW_MERGE(EPT, MIN(EPT, 2*EPA+1-J), DZ, I, J, DZ, K, L)
       ELSE IF (T .EQ. 0) THEN
          L = EPT - 1
          DO J = 0, L
             DZ(K+J) = DZ(I+J)
          END DO
          L = 0
       ELSE ! nothing to do
          L = 0
       END IF
       !$OMP END PARALLEL
       TE = TE + L

       IF (TE .EQ. 0) EXIT
    END DO

#ifndef USE_FAST
    INFO = GET_SYS_TIME() - INFO
#endif
  END SUBROUTINE AW_SRT1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE AW_SRT2(NN, NM, DZ, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NN, NM
    TYPE(AW), INTENT(INOUT) :: DZ(NM)
    INTEGER, INTENT(OUT) :: INFO
    
    IF (NN .LT. 0) THEN
       INFO = -1
    ELSE IF (NM .LT. NN) THEN
       INFO = -2
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN
    IF (NN .EQ. 0) RETURN

#ifndef USE_FAST
    INFO = GET_SYS_TIME()
#endif
    CALL AW_SORT(NN, DZ, DZ(NN+1))
#ifndef USE_FAST
    INFO = GET_SYS_TIME() - INFO
#endif
  END SUBROUTINE AW_SRT2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE AW_NCP1(N, NN, NM, DZ, N_2, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, NN, NM, N_2
    TYPE(AW), INTENT(INOUT) :: DZ(NM)
    INTEGER, INTENT(OUT) :: SL, STEP(N_2), INFO

    INTEGER :: J, BP, BQ
    INTEGER :: I, K, AP, AQ
    !$ SAVE :: I, K, AP, AQ

    SL = 0

    IF (N .LT. 0) THEN
       INFO = -1
    ELSE IF (NN .LT. 0) THEN
       INFO = -2
    ELSE IF (NM .LT. NN) THEN
       INFO = -3
    ELSE IF (N_2 .LT. 0) THEN
       INFO = -5
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

#ifndef USE_FAST
    INFO = GET_SYS_TIME()
#endif
    IF (N .EQ. 0) GOTO 1
    IF (NN .EQ. 0) GOTO 1
    IF (N_2 .EQ. 0) GOTO 1

    I = 1
    DO WHILE (SL .LT. N_2)
       SL = SL + 1
       STEP(SL) = I
       IF (SL .GE. N_2) EXIT

       AP = DZ(I)%P
       AQ = DZ(I)%Q
       K = NN + 1

       !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J,BP,BQ) SHARED(I,NN,DZ,AP,AQ) REDUCTION(MIN:K)
       DO J = I+1, NN
          IF (DZ(J)%W .EQ. DZ(J)%W) THEN
             BP = DZ(J)%P
             BQ = DZ(J)%Q
             IF ((AP .EQ. BP) .OR. (AP .EQ. BQ) .OR. (AQ .EQ. BP) .OR. (AQ .EQ. BQ)) THEN
                DZ(J)%W = QUIET_NAN(J)
             ELSE ! not colliding
                K = MIN(K, J)
             END IF
          END IF
       END DO
       !$OMP END PARALLEL DO

       IF (K .GT. NN) EXIT
       I = K
    END DO

1   CONTINUE
#ifndef USE_FAST
    INFO = GET_SYS_TIME() - INFO
#endif
  END SUBROUTINE AW_NCP1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE AW_NCP2(N, NN, NM, DZ, N_2, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, NN, NM, N_2
    TYPE(AW), INTENT(INOUT) :: DZ(NM)
    INTEGER, INTENT(OUT) :: SL, STEP(N_2), INFO

    LOGICAL :: R(N)
    INTEGER :: I, J

    SL = 0

    IF (N .LT. 0) THEN
       INFO = -1
    ELSE IF (NN .LT. 0) THEN
       INFO = -2
    ELSE IF (NM .LT. NN) THEN
       INFO = -3
    ELSE IF (N_2 .LT. 0) THEN
       INFO = -5
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

#ifndef USE_FAST
    INFO = GET_SYS_TIME()
#endif
    IF (N .EQ. 0) GOTO 2
    IF (NN .EQ. 0) GOTO 2
    IF (N_2 .EQ. 0) GOTO 2

    R = .FALSE.
    J = 1

    DO I = 1, N_2
       STEP(I) = J
       R(DZ(J)%P) = .TRUE.
       R(DZ(J)%Q) = .TRUE.
       J = J + 1
       DO WHILE (J .LE. NN)
          IF ((DZ(J)%W .EQ. DZ(J)%W) .AND. (.NOT. R(DZ(J)%P)) .AND. (.NOT. R(DZ(J)%Q))) EXIT
          J = J + 1
       END DO
       SL = I
       IF (J .GT. NN) EXIT
    END DO

2   CONTINUE
#ifndef USE_FAST
    INFO = GET_SYS_TIME() - INFO
#endif
  END SUBROUTINE AW_NCP2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE AW_NCPT(W, J, N, NN, NM, DZ, N_2, SL, STEP, K)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(INOUT) :: W
    INTEGER, INTENT(IN) :: J, N, NN, NM, N_2
    TYPE(AW), INTENT(INOUT) :: DZ(NM)
    INTEGER, INTENT(INOUT) :: SL, STEP(N_2), K

    INTEGER :: STP(N_2)
    REAL(KIND=DWP) :: MYW, MYWP, MYWN
    INTEGER :: I, L, MYSL

    CALL AW_NCP2(N, NN, NM, DZ, N_2, MYSL, STP, I)
    IF (I .LT. 0) RETURN

    IF (MYSL .GE. 1) THEN
       L = 0
       DO I = MYSL, 1, -1
          IF (.NOT. (DZ(STP(I))%W .LT. D_ZERO)) THEN
             L = I
             EXIT
          END IF
       END DO

       MYWP = D_ZERO
       DO I = L, 1, -1
          MYWP = MYWP + DZ(STP(I))%W
       END DO

       MYWN = D_ZERO
       DO I = L+1, MYSL
          MYWN = MYWN + DZ(STP(I))%W
       END DO

       MYW = MYWP + MYWN
       L = J + 1
       !$OMP CRITICAL
       IF ((.NOT. (MYW .LE. W)) .OR. ((MYW .EQ. W) .AND. (L .LT. K))) THEN
          SL = MYSL
          DO I = 1, SL
             STEP(I) = STP(I) + J
          END DO
          W = MYW
          K = L
       END IF
       !$OMP END CRITICAL
    END IF
  END SUBROUTINE AW_NCPT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE AW_NCP3(N, NN, NM, DZ, N_2, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, NN, NM, N_2
    TYPE(AW), INTENT(INOUT) :: DZ(NM)
    INTEGER, INTENT(OUT) :: SL, STEP(N_2), INFO

    INTEGER :: I, J
    REAL(KIND=DWP) :: W
    INTEGER :: K
    !$ SAVE :: W, K

    SL = 0

    IF (N .LT. 0) THEN
       INFO = -1
    ELSE IF (NN .LT. 0) THEN
       INFO = -2
    ELSE IF (NM .LT. NN) THEN
       INFO = -3
    ELSE IF (N_2 .LT. 0) THEN
       INFO = -5
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

#ifndef USE_FAST
    INFO = GET_SYS_TIME()
#endif
    IF (N .EQ. 0) GOTO 3
    IF (NN .EQ. 0) GOTO 3
    IF (N_2 .EQ. 0) GOTO 3

    W = QUIET_NAN(N_2)
    K = NN + 1
#ifdef USE_INTEL
    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) SHARED(W,K,N,NN,NM,DZ,N_2,SL,STEP) SCHEDULE(NONMONOTONIC:DYNAMIC,1)
#else
    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,J) SHARED(W,K,N,NN,NM,DZ,N_2,SL,STEP) SCHEDULE(DYNAMIC,1)
#endif
    DO I = 1, NN
       J = I - 1
       CALL AW_NCPT(W, J, N, NN - J, NM - J, DZ(I), N_2, SL, STEP, K)
    END DO
    !$OMP END PARALLEL DO

3   CONTINUE
#ifndef USE_FAST
    INFO = GET_SYS_TIME() - INFO
#endif
  END SUBROUTINE AW_NCP3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE atypes
