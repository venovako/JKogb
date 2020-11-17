MODULE BINIO
  USE PARAMS
  IMPLICIT NONE

  INTEGER, PARAMETER, PRIVATE :: IOMSGLEN = 66

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE BIO_OPEN(U, FN, ACT, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: U
    CHARACTER(LEN=*), INTENT(IN) :: FN, ACT
    INTEGER, INTENT(OUT) :: INFO

    CHARACTER(LEN=IOMSGLEN) :: MSG
    CHARACTER(LEN=9) :: FACT
    CHARACTER(LEN=7) :: STAT

    IF (U .LT. -1) THEN
       INFO = -1
    ELSE IF (LEN_TRIM(FN) .LE. 0) THEN
       INFO = -2
    ELSE IF (LEN_TRIM(ACT) .NE. 2) THEN
       INFO = -3
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    SELECT CASE (ACT(1:1))
    CASE ('R', 'r')
       STAT = 'OLD'
       SELECT CASE (ACT(2:2))
       CASE ('O', 'o')
          FACT = 'READ'
       CASE ('W', 'w')
          FACT = 'READWRITE'
       CASE DEFAULT
          INFO = -3
       END SELECT
    CASE ('W', 'w')
       STAT = 'REPLACE'
       SELECT CASE (ACT(2:2))
       CASE ('O', 'o')
          FACT = 'WRITE'
       CASE ('R', 'r')
          FACT = 'READWRITE'
       CASE DEFAULT
          INFO = -3
       END SELECT
    CASE DEFAULT
       INFO = -3
    END SELECT
    IF (INFO .NE. 0) RETURN

    IF (U .EQ. -1) THEN
       OPEN(NEWUNIT=U, IOSTAT=INFO, IOMSG=MSG, ERR=1, FILE=TRIM(FN), STATUS=TRIM(STAT), ACTION=TRIM(FACT),&
            ACCESS='STREAM', FORM='UNFORMATTED')
    ELSE
       OPEN(UNIT=U, IOSTAT=INFO, IOMSG=MSG, ERR=1, FILE=TRIM(FN), STATUS=TRIM(STAT), ACTION=TRIM(FACT),&
            ACCESS='STREAM', FORM='UNFORMATTED')
    END IF
    RETURN
1   WRITE (ERROR_UNIT,'(2A)') 'BIO_OPEN    :', TRIM(MSG)
  END SUBROUTINE BIO_OPEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE BIO_CLOSE(U, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: U
    INTEGER, INTENT(OUT) :: INFO

    CHARACTER(LEN=IOMSGLEN) :: MSG

    INFO = 0
    CLOSE(UNIT=U, IOSTAT=INFO, IOMSG=MSG, ERR=2)
    U = -1
    RETURN
2   WRITE (ERROR_UNIT,'(2A)') 'BIO_CLOSE   :', TRIM(MSG)
  END SUBROUTINE BIO_CLOSE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE BIO_READ_I1(U, M, J, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: U, M
    INTEGER, INTENT(OUT) :: J(M), INFO

    CHARACTER(LEN=IOMSGLEN) :: MSG
    INTEGER(KIND=DWP) :: K
    INTEGER :: I

    IF (M .LT. 0) THEN
       INFO = -2
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    DO I = 1, M
       READ (UNIT=U, IOSTAT=INFO, IOMSG=MSG, ERR=3) K
       J(I) = INT(K)
    END DO
    RETURN
3   WRITE (ERROR_UNIT,'(2A)') 'BIO_READ_I1 :', TRIM(MSG)
  END SUBROUTINE BIO_READ_I1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE BIO_WRITE_D1(U, N, A, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: U, N
    REAL(KIND=DWP), INTENT(IN) :: A(N)
    INTEGER, INTENT(OUT) :: INFO

    CHARACTER(LEN=IOMSGLEN) :: MSG
    INTEGER :: I

    IF (N .LT. 0) THEN
       INFO = -2
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    DO I = 1, N
       WRITE (UNIT=U, IOSTAT=INFO, IOMSG=MSG, ERR=4) A(I)
    END DO
    RETURN
4   WRITE (ERROR_UNIT,'(2A)') 'BIO_WRITE_D1:', TRIM(MSG)
  END SUBROUTINE BIO_WRITE_D1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE BIO_READ_D2(U, M, N, A, LDA, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: U, M, N, LDA
    REAL(KIND=DWP), INTENT(OUT) :: A(LDA,N)
    INTEGER, INTENT(OUT) :: INFO

    CHARACTER(LEN=IOMSGLEN) :: MSG
    INTEGER :: I, J

    IF (M .LT. 0) THEN
       INFO = -2
    ELSE IF (N .LT. 0) THEN
       INFO = -3
    ELSE IF (LDA .LT. M) THEN
       INFO = -5
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    DO J = 1, N
       DO I = 1, M
          READ (UNIT=U, IOSTAT=INFO, IOMSG=MSG, ERR=5) A(I,J)
       END DO
    END DO
    RETURN
5   WRITE (ERROR_UNIT,'(2A)') 'BIO_READ_D2 :', TRIM(MSG)
  END SUBROUTINE BIO_READ_D2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE BIO_WRITE_D2(U, M, N, A, LDA, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: U, M, N, LDA
    REAL(KIND=DWP), INTENT(IN) :: A(LDA,N)
    INTEGER, INTENT(OUT) :: INFO

    CHARACTER(LEN=IOMSGLEN) :: MSG
    INTEGER :: I, J

    IF (M .LT. 0) THEN
       INFO = -2
    ELSE IF (N .LT. 0) THEN
       INFO = -3
    ELSE IF (LDA .LT. M) THEN
       INFO = -5
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    DO J = 1, N
       DO I = 1, M
          WRITE (UNIT=U, IOSTAT=INFO, IOMSG=MSG, ERR=6) A(I,J)
       END DO
    END DO
    RETURN
6   WRITE (ERROR_UNIT,'(2A)') 'BIO_WRITE_D2:', TRIM(MSG)
  END SUBROUTINE BIO_WRITE_D2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE BIO_READ_Z2(U, M, N, A, LDA, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: U, M, N, LDA
    COMPLEX(KIND=DWP), INTENT(OUT) :: A(LDA,N)
    INTEGER, INTENT(OUT) :: INFO

    CHARACTER(LEN=IOMSGLEN) :: MSG
    INTEGER :: I, J

    IF (M .LT. 0) THEN
       INFO = -2
    ELSE IF (N .LT. 0) THEN
       INFO = -3
    ELSE IF (LDA .LT. M) THEN
       INFO = -5
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    DO J = 1, N
       DO I = 1, M
          READ (UNIT=U, IOSTAT=INFO, IOMSG=MSG, ERR=7) A(I,J)
       END DO
    END DO
    RETURN
7   WRITE (ERROR_UNIT,'(2A)') 'BIO_READ_Z2 :', TRIM(MSG)
  END SUBROUTINE BIO_READ_Z2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE BIO_WRITE_Z2(U, M, N, A, LDA, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: U, M, N, LDA
    COMPLEX(KIND=DWP), INTENT(IN) :: A(LDA,N)
    INTEGER, INTENT(OUT) :: INFO

    CHARACTER(LEN=IOMSGLEN) :: MSG
    INTEGER :: I, J

    IF (M .LT. 0) THEN
       INFO = -2
    ELSE IF (N .LT. 0) THEN
       INFO = -3
    ELSE IF (LDA .LT. M) THEN
       INFO = -5
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    DO J = 1, N
       DO I = 1, M
          WRITE (UNIT=U, IOSTAT=INFO, IOMSG=MSG, ERR=8) A(I,J)
       END DO
    END DO
    RETURN
8   WRITE (ERROR_UNIT,'(2A)') 'BIO_WRITE_Z2:', TRIM(MSG)
  END SUBROUTINE BIO_WRITE_Z2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE OPEN_YJ_RO(FN, FD, INFO)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: FN
    INTEGER, INTENT(OUT) :: FD(2), INFO

    FD = -1

    IF (LEN_TRIM(FN) .LE. 0) THEN
       INFO = -1
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    CALL BIO_OPEN(FD(1), (TRIM(FN)//'.Y'), 'RO', INFO)
    IF (INFO .NE. 0) THEN
       INFO = 1
       RETURN
    END IF

    CALL BIO_OPEN(FD(2), (TRIM(FN)//'.J'), 'RO', INFO)
    IF (INFO .NE. 0) THEN
       INFO = 2
       RETURN
    END IF
  END SUBROUTINE OPEN_YJ_RO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE OPEN_UZS_WO(FN, FD, INFO)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: FN
    INTEGER, INTENT(OUT) :: FD(3), INFO

    FD = -1

    IF (LEN_TRIM(FN) .LE. 0) THEN
       INFO = -1
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    CALL BIO_OPEN(FD(1), (TRIM(FN)//'.YU'), 'WO', INFO)
    IF (INFO .NE. 0) THEN
       INFO = 1
       RETURN
    END IF

    CALL BIO_OPEN(FD(2), (TRIM(FN)//'.ZZ'), 'WO', INFO)
    IF (INFO .NE. 0) THEN
       INFO = 2
       RETURN
    END IF

    CALL BIO_OPEN(FD(3), (TRIM(FN)//'.SY'), 'WO', INFO)
    IF (INFO .NE. 0) THEN
       INFO = 3
       RETURN
    END IF
  END SUBROUTINE OPEN_UZS_WO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE BCLOSEN(N, FD, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(INOUT) :: FD(N)
    INTEGER, INTENT(OUT) :: INFO

    INTEGER :: I

    IF (N .LT. 0) THEN
       INFO = -1
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    DO I = N, 1, -1
       CALL BIO_CLOSE(FD(I), INFO)
       IF (INFO .NE. 0) THEN
          INFO = I
          RETURN
       END IF
    END DO
  END SUBROUTINE BCLOSEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE BINIO
