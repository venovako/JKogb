!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RECURSIVE SUBROUTINE AW_NCP2(NT, N, NN, NM, DZ, N_2, SL, STEP, INFO)
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

    !DIR$ VECTOR ALWAYS
    DO I = 1, N
       P(I) = .FALSE.
    END DO
    !DIR$ VECTOR ALWAYS
    DO I = 1, N
       Q(I) = .FALSE.
    END DO

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

2   INFO = GET_SYS_US() - INFO
  END SUBROUTINE AW_NCP2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RECURSIVE SUBROUTINE AW_NCPT(W, N, NN, NM, DZ, N_2, SL, STEP)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(INOUT) :: W
    INTEGER, INTENT(IN) :: N, NN, NM, N_2
    TYPE(AW), INTENT(INOUT) :: DZ(NM)
    INTEGER, INTENT(INOUT) :: SL, STEP(N_2)

    INTEGER, TARGET :: STP(N_2)
    REAL(KIND=DWP) :: MYW
    INTEGER :: I, MYSL

    CALL AW_NCP2(1, N, NN, NM, DZ, N_2, MYSL, STP, I)
    IF (I .LT. 0) RETURN

    MYW = D_ZERO
    DO I = SL, 1, -1
       MYW = MYW + DZ(STP(I))%W
    END DO

    !$OMP CRITICAL
    IF (.NOT. (MYW .LE. W)) THEN
       SL = MYSL
       DO I = 1, SL
          STEP(I) = STP(I)
       END DO
       W = MYW
    END IF
    !$OMP END CRITICAL
  END SUBROUTINE AW_NCPT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RECURSIVE SUBROUTINE AW_NCP3(NT, N, NN, NM, DZ, N_2, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, N, NN, NM, N_2
    TYPE(AW), INTENT(INOUT) :: DZ(NM)
    INTEGER, INTENT(OUT) :: SL, STEP(N_2), INFO

    REAL(KIND=DWP) :: W
    INTEGER :: I

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

    W = QUIET_NAN(N_2)
    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I) SHARED(W,N,NN,NM,DZ,N_2,SL,STEP)
    DO I = 1, NN
       CALL AW_NCPT(W, N, NN - I + 1, NM - I + 1, DZ(I), N_2, SL, STEP)
    END DO
    !$OMP END PARALLEL DO

3   INFO = GET_SYS_US() - INFO
  END SUBROUTINE AW_NCP3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
