!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RECURSIVE SUBROUTINE AW_NCP4(NT, N, NN, NM, DZ, N_2, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, N, NN, NM, N_2
    TYPE(AW), INTENT(INOUT) :: DZ(NM)
    INTEGER, INTENT(OUT) :: SL, STEP(N_2), INFO

    LOGICAL, TARGET :: P(N), Q(N)
    INTEGER :: I, J

    SL = 0
#ifndef NDEBUG
    DO I = 1, N_2
       STEP(I) = 0
    END DO
#endif

    IF (NT .NE. 1) THEN
       INFO = -1
    ELSE IF (N .LT. 0) THEN
       INFO = -2
    ELSE IF (NN .LT. 0) THEN
       INFO = -3
    ELSE IF (NM .LT. NN) THEN
       INFO = -4
    ELSE IF (N_2 .LT. 0) THEN
       INFO = -6
    ELSE IF (N_2 .GT. MIN(N,NN)) THEN
       INFO = -6
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    INFO = GET_SYS_US()
    IF (N .EQ. 0) GOTO 4
    IF (NN .EQ. 0) GOTO 4
    IF (N_2 .EQ. 0) GOTO 4

    P = .FALSE.
    Q = .FALSE.
    J = 1

    DO I = 1, N_2
       STEP(I) = J
#ifndef NDEBUG
       IF (P(DZ(J)%P)) THEN
          INFO = -5
          RETURN
       ELSE
#endif
          P(DZ(J)%P) = .TRUE.
#ifndef NDEBUG
       END IF
       IF (Q(DZ(J)%Q) THEN
          INFO = -5
          RETURN
       ELSE
#endif
          Q(DZ(J)%Q) = .TRUE.
#ifndef NDEBUG
       END IF
#endif
       J = J + 1
       DO WHILE (J .LE. NN)
          IF ((DZ(J)%W .EQ. DZ(J)%W) .AND. (.NOT. P(DZ(J)%P)) .AND. (.NOT. Q(DZ(J)%Q))) EXIT
          J = J + 1
       END DO
       SL = I
       IF (J .GT. NN) EXIT
    END DO

4   INFO = GET_SYS_US() - INFO
  END SUBROUTINE AW_NCP4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RECURSIVE SUBROUTINE AW_NCP5(NT, N, NN, NM, DZ, N_2, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, N, NN, NM, N_2
    TYPE(AW), INTENT(INOUT) :: DZ(NM)
    INTEGER, INTENT(OUT) :: SL, STEP(N_2), INFO

    INTEGER, TARGET :: STP(N_2,NT)
    INTEGER :: I, J, K, L, T
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
    ELSE IF (N_2 .GT. MIN(N,NN)) THEN
       INFO = -6
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    INFO = GET_SYS_US()
    IF (N .EQ. 0) GOTO 5
    IF (NN .EQ. 0) GOTO 5
    IF (N_2 .EQ. 0) GOTO 5

    W1 = QUIET_NAN(N_2)
    J = 0
    !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(I,K,L,T,WL) SHARED(N,NN,NM,DZ,N_2,SL,STEP,STP,J,W1)
    DO K = 1, NN
       T = INT(OMP_GET_THREAD_NUM()) + 1
       CALL AW_NCP2(1, N, NN, NM, DZ, N_2, L, STP(:,T), I)
#ifndef NDEBUG
       IF (I .LT. 0) STOP
#endif
       WL = D_ZERO
       DO I = L, 1, -1
          WL = WL + DZ(STP(I,T))%W
       END DO
       !$OMP CRITICAL
       IF (.NOT. (WL .LE. W1)) THEN
          SL = L
          DO I = 1, SL
             STEP(I) = STP(I,T)
          END DO
          J = K
          W1 = WL
       END IF
       !$OMP END CRITICAL
    END DO
    !$OMP END PARALLEL DO

5   INFO = GET_SYS_US() - INFO
  END SUBROUTINE AW_NCP5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
