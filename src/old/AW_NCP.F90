!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE AW_NCP(NN, DZ, N_2, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NN, N_2
    TYPE(AW), INTENT(IN) :: DZ(NN)
    INTEGER, INTENT(OUT) :: SL, STEP(N_2), INFO

    INTEGER :: I, J, K, AP, AQ, BP, BQ
    LOGICAL :: C

    IF (NN .LT. 0) THEN
       INFO = -1
    ELSE IF (N_2 .LT. 0) THEN
       INFO = -3
    ELSE IF (N_2 .GT. NN) THEN
       INFO = -3
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    INFO = GET_THREAD_NS()
    SL = 0
    IF (NN .EQ. 0) GOTO 1
    IF (N_2 .EQ. 0) GOTO 1

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

    !DIR$ VECTOR ALWAYS
    DO I = SL+1, N_2
       STEP(I) = 0
    END DO

1   INFO = GET_THREAD_NS() - INFO
  END SUBROUTINE AW_NCP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
