MODULE JSTEP
  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! column-cyclic
  PURE SUBROUTINE TRU1(N, J, NN, P, Q, PN, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, J(N), NN
    INTEGER, INTENT(OUT) :: P(NN), Q(NN), PN(2), INFO

    INTEGER :: IP, IQ, IT, IH

    PN = 0

    IF (N .LT. 0) THEN
       INFO = -1
    ELSE IF (NN .LT. 0) THEN
       INFO = -3
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN
    IF (N .LE. 1) RETURN

    DO IT = 1, N
       SELECT CASE (J(IT))
       CASE (-1)
          PN(2) = PN(2) + 1
       CASE (1)
          PN(1) = PN(1) + 1
       CASE DEFAULT
          INFO = -2
          RETURN
       END SELECT
    END DO

#ifndef NDEBUG
    DO IT = 1, NN
       P(IT) = 0
    END DO

    DO IT = 1, NN
       Q(IT) = 0
    END DO
#endif

    IT = 1
    IH = NN
    DO IQ = 2, N
       DO IP = 1, IQ-1
          IF (IT .GT. IH) RETURN
          IF (J(IP) .EQ. J(IQ)) THEN
             P(IT) = IP
             Q(IT) = IQ
             IT = IT + 1
          ELSE ! hyp
             P(IH) = IP
             Q(IH) = IQ
             IH = IH - 1
          END IF
       END DO
    END DO

#ifndef NDEBUG
    DO IH = 1, NN
       IF (P(IH) .EQ. 0) THEN
          INFO = -4
          RETURN
       END IF
    END DO

    DO IH = 1, NN
       IF (Q(IH) .EQ. 0) THEN
          INFO = -5
          RETURN
       END IF
    END DO
#endif
    INFO = IT - 1
  END SUBROUTINE TRU1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! row-cyclic
  PURE SUBROUTINE TRU2(N, J, NN, P, Q, PN, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, J(N), NN
    INTEGER, INTENT(OUT) :: P(NN), Q(NN), PN(2), INFO

    INTEGER :: IP, IQ, IT, IH

    PN = 0

    IF (N .LT. 0) THEN
       INFO = -1
    ELSE IF (NN .LT. 0) THEN
       INFO = -3
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN
    IF (N .LE. 1) RETURN

    DO IT = 1, N
       SELECT CASE (J(IT))
       CASE (-1)
          PN(2) = PN(2) + 1
       CASE (1)
          PN(1) = PN(1) + 1
       CASE DEFAULT
          INFO = -2
          RETURN
       END SELECT
    END DO

#ifndef NDEBUG
    DO IT = 1, NN
       P(IT) = 0
    END DO

    DO IT = 1, NN
       Q(IT) = 0
    END DO
#endif

    IT = 1
    IH = NN
    DO IP = 1, N-1
       DO IQ = IP+1, N
          IF (IT .GT. IH) RETURN
          IF (J(IP) .EQ. J(IQ)) THEN
             P(IT) = IP
             Q(IT) = IQ
             IT = IT + 1
          ELSE ! hyp
             P(IH) = IP
             Q(IH) = IQ
             IH = IH - 1
          END IF
       END DO
    END DO

#ifndef NDEBUG
    DO IH = 1, NN
       IF (P(IH) .EQ. 0) THEN
          INFO = -4
          RETURN
       END IF
    END DO

    DO IH = 1, NN
       IF (Q(IH) .EQ. 0) THEN
          INFO = -5
          RETURN
       END IF
    END DO
#endif
    INFO = IT - 1
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
