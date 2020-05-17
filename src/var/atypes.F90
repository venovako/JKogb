!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RECURSIVE SUBROUTINE AW_NCP2(NT, N, NN, NM, DZ, N_2, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, N, NN, NM, N_2
    TYPE(AW), INTENT(INOUT) :: DZ(NM)
    INTEGER, INTENT(OUT) :: SL, STEP(N_2), INFO

    INTEGER :: I, J, K, AP, AQ, BP, BQ
    LOGICAL :: C

    SL = 0
#ifndef NDEBUG
    DO I = 1, N_2
       STEP(I) = 0
    END DO
#endif

    IF (NT .LE. 0) THEN
       INFO = -1
    ELSE IF (N .LT. 0) THEN
       INFO = -2
    ELSE IF (NN .LT. 0) THEN
       INFO = -3
    ELSE IF (NM .LT. NN) THEN
       INFO = -4
    ELSE IF (N_2 .LT. 0) THEN
       INFO = -6
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    INFO = GET_SYS_US()
    IF (N .EQ. 0) GOTO 2
    IF (NN .EQ. 0) GOTO 2
    IF (N_2 .EQ. 0) GOTO 2

    J = 1
    DO I = 1, N_2
       STEP(I) = J
       J = J + 1
       DO WHILE (J .LE. NN)
          C = .NOT. (DZ(J)%W .EQ. DZ(J)%W)
          IF (C) THEN
             J = J + 1
             CYCLE
          END IF
          AP = DZ(J)%P
          AQ = DZ(J)%Q
          DO K = I, 1, -1
             BP = DZ(STEP(K))%P
             BQ = DZ(STEP(K))%Q
             IF ((AP .EQ. BP) .OR. (AP .EQ. BQ) .OR. (AQ .EQ. BP) .OR. (AQ .EQ. BQ)) THEN
                C = .TRUE.
                EXIT
             END IF
          END DO
          IF (C) THEN
             J = J + 1
          ELSE
             EXIT
          END IF
       END DO
       SL = I
       IF (J .GT. NN) EXIT
    END DO

2   INFO = GET_SYS_US() - INFO
  END SUBROUTINE AW_NCP2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RECURSIVE SUBROUTINE AW_NCP3(NT, NN, N, NM, DZ, N_2, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, N, NN, NM, N_2
    TYPE(AW), INTENT(INOUT) :: DZ(NM)
    INTEGER, INTENT(OUT) :: SL, STEP(N_2), INFO

    INTEGER, TARGET :: STP(N_2)
    INTEGER :: I, J, K, L, M, AP, AQ, BP, BQ
    LOGICAL :: C
    REAL(KIND=DWP) :: W1, WL

    SL = 0
#ifndef NDEBUG
    DO I = 1, N_2
       STEP(I) = 0
    END DO
#endif

    IF (NT .LE. 0) THEN
       INFO = -1
    ELSE IF (N .LT. 0) THEN
       INFO = -2
    ELSE IF (NN .LT. 0) THEN
       INFO = -3
    ELSE IF (NM .LT. NN) THEN
       INFO = -4
    ELSE IF (N_2 .LT. 0) THEN
       INFO = -6
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    INFO = GET_SYS_US()
    IF (N .EQ. 0) GOTO 3
    IF (NN .EQ. 0) GOTO 3
    IF (N_2 .EQ. 0) GOTO 3

    J = 1
    DO I = 1, N_2
       STEP(I) = J
       J = J + 1
       DO WHILE (J .LE. NN)
          C = .NOT. (DZ(J)%W .EQ. DZ(J)%W)
          IF (C) THEN
             J = J + 1
             CYCLE
          END IF
          AP = DZ(J)%P
          AQ = DZ(J)%Q
          DO K = I, 1, -1
             BP = DZ(STEP(K))%P
             BQ = DZ(STEP(K))%Q
             IF ((AP .EQ. BP) .OR. (AP .EQ. BQ) .OR. (AQ .EQ. BP) .OR. (AQ .EQ. BQ)) THEN
                C = .TRUE.
                EXIT
             END IF
          END DO
          IF (C) THEN
             J = J + 1
          ELSE
             EXIT
          END IF
       END DO
       SL = I
       IF (J .GT. NN) EXIT
    END DO
    IF (SL .EQ. 0) RETURN
    W1 = D_ZERO
    DO I = SL, 1, -1
       W1 = W1 + DZ(STEP(I))%W
    END DO

#ifdef OLD_OMP
    !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(I,J,K,L,M,AP,AQ,BP,BQ,C,WL,STP) SHARED(NN,N_2,SL,DZ,STEP,W1) &
    !$OMP& SCHEDULE(DYNAMIC,1)
#else
    !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(I,J,K,L,M,AP,AQ,BP,BQ,C,WL,STP) SHARED(NN,N_2,SL,DZ,STEP,W1) &
    !$OMP& SCHEDULE(NONMONOTONIC:DYNAMIC,1)
#endif
    DO L = 2, NN
       M = 0
       J = L
       DO I = 1, N_2
          STP(I) = J
          J = J + 1
          DO WHILE (J .LE. NN)
             C = .NOT. (DZ(J)%W .EQ. DZ(J)%W)
             IF (C) THEN
                J = J + 1
                CYCLE
             END IF
             AP = DZ(J)%P
             AQ = DZ(J)%Q
             DO K = I, 1, -1
                BP = DZ(STP(K))%P
                BQ = DZ(STP(K))%Q
                IF ((AP .EQ. BP) .OR. (AP .EQ. BQ) .OR. (AQ .EQ. BP) .OR. (AQ .EQ. BQ)) THEN
                   C = .TRUE.
                   EXIT
                END IF
             END DO
             IF (C) THEN
                J = J + 1
             ELSE
                EXIT
             END IF
          END DO
          M = I
          IF (J .GT. NN) EXIT
       END DO
       WL = D_ZERO
       DO I = M, 1, -1
          WL = WL + DZ(STP(I))%W
       END DO
       !$OMP CRITICAL
       IF (.NOT. (WL .LE. W1)) THEN
          SL = M
          DO I = 1, SL
             STEP(I) = STP(I)
          END DO
          W1 = WL
       END IF
       !$OMP END CRITICAL
    END DO
    !$OMP END PARALLEL DO

3   INFO = GET_SYS_US() - INFO
  END SUBROUTINE AW_NCP3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
