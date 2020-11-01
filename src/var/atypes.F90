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

  RECURSIVE SUBROUTINE AW_NCPT(W, J, N, NN, NM, DZ, N_2, SL, STEP, K)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(INOUT) :: W
    INTEGER, INTENT(IN) :: J, N, NN, NM, N_2
    TYPE(AW), INTENT(INOUT) :: DZ(NM)
    INTEGER, INTENT(INOUT) :: SL, STEP(N_2), K

    INTEGER, TARGET :: STP(N_2)
    REAL(KIND=DWP) :: MYW, MYWP, MYWN
    INTEGER :: I, L, MYSL

    CALL AW_NCP2(1, N, NN, NM, DZ, N_2, MYSL, STP, I)
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

  RECURSIVE SUBROUTINE AW_NCP3(NT, N, NN, NM, DZ, N_2, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, N, NN, NM, N_2
    TYPE(AW), INTENT(INOUT) :: DZ(NM)
    INTEGER, INTENT(OUT) :: SL, STEP(N_2), INFO

    REAL(KIND=DWP) :: W
    INTEGER :: I, J, K

    SL = 0

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

3   INFO = GET_SYS_US() - INFO
  END SUBROUTINE AW_NCP3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
