MODULE JSTEP
  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! column-cyclic
  PURE INTEGER FUNCTION PQI1(N, P, Q)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, P, Q

    IF (N .LT. 2) THEN
       PQI1 = -1
    ELSE IF ((P .LT. 1) .OR. (P .GT. N)) THEN
       PQI1 = -2
    ELSE IF ((Q .LT. 2) .OR. (Q .GT. N)) THEN
       PQI1 = -3
    ELSE IF (P .GE. Q) THEN
       PQI1 = 0
    ELSE ! all OK
       PQI1 = ((Q - 2) * (Q - 1)) / 2 + P
    END IF
  END FUNCTION PQI1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! row-cyclic
  PURE INTEGER FUNCTION PQI2(N, P, Q)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, P, Q

    IF (N .LT. 2) THEN
       PQI2 = -1
    ELSE IF ((P .LT. 1) .OR. (P .GT. N)) THEN
       PQI2 = -2
    ELSE IF ((Q .LT. 2) .OR. (Q .GT. N)) THEN
       PQI2 = -3
    ELSE IF (P .GE. Q) THEN
       PQI2 = 0
    ELSE ! all OK
       PQI2 = (2 * N - P) * (P - 1) + (Q - P)
    END IF
  END FUNCTION PQI2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! column-cyclic
  PURE SUBROUTINE TRU1(N, P, Q, NN, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, NN
    INTEGER, INTENT(OUT) :: P(NN), Q(NN), INFO

    INTEGER :: IP, IQ, I

    IF (N .LT. 0) THEN
       INFO = -1
    ELSE IF (NN .LT. 0) THEN
       INFO = -4
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN
    IF (N .LE. 1) RETURN

    INFO = (N * (N - 1)) / 2

    I = 1
    DO IQ = 2, N
       DO IP = 1, IQ-1
          IF (I .GT. NN) RETURN
          P(I) = IP
          Q(I) = IQ
          I = I + 1
       END DO
    END DO
  END SUBROUTINE TRU1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! row-cyclic
  PURE SUBROUTINE TRU2(N, P, Q, NN, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, NN
    INTEGER, INTENT(OUT) :: P(NN), Q(NN), INFO

    INTEGER :: IP, IQ, I

    IF (N .LT. 0) THEN
       INFO = -1
    ELSE IF (NN .LT. 0) THEN
       INFO = -4
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN
    IF (N .LE. 1) RETURN

    INFO = (N * (N - 1)) / 2

    I = 1
    DO IP = 1, N-1
       DO IQ = IP+1, N
          IF (I .GT. NN) RETURN
          P(I) = IP
          Q(I) = IQ
          I = I + 1
       END DO
    END DO
  END SUBROUTINE TRU2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE INTEGER FUNCTION JSTEP_LEN(N, N_2)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, N_2

    IF (N_2 .EQ. 0) THEN
       JSTEP_LEN = MAX((N / 2), 0)
    ELSE
       JSTEP_LEN = MIN(MAX(N_2, 0), (N / 2))
    END IF
  END FUNCTION JSTEP_LEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE JSTEP
