!>@brief \b BFOPEN opens a file for binary I/O and returns its newly assigned unit number.
!!
!!@param FN [IN]; a file name.
!!@param ACT [IN]; one of RO, RW (or WR), or WO case-insensitive strings for Read-Only, Read-Write, or Write-Only access mode, respectively.
!!@param U [OUT]; the file's newly assigned unit number.
!!@param INFO [OUT]; zero on success, \f$-i\f$ if the \f$i\f$th argument had an illegal value, or a positive I/O error code.
SUBROUTINE BFOPEN(FN, ACT, U, INFO)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: FN
  CHARACTER(LEN=2), INTENT(IN) :: ACT
  INTEGER, INTENT(OUT) :: U, INFO

  CHARACTER(LEN=9) :: FACT
  CHARACTER(LEN=7) :: STAT

  IF (LEN_TRIM(FN) .LE. 0) THEN
     INFO = -1
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
        INFO = -2
     END SELECT
  CASE ('W', 'w')
     STAT = 'REPLACE'
     SELECT CASE (ACT(2:2))
     CASE ('O', 'o')
        FACT = 'WRITE'
     CASE ('R', 'r')
        FACT = 'READWRITE'
     CASE DEFAULT
        INFO = -2
     END SELECT
  CASE DEFAULT
     INFO = -2
  END SELECT
  IF (INFO .NE. 0) RETURN

  OPEN(NEWUNIT=U, IOSTAT=INFO, FILE=TRIM(FN), STATUS=TRIM(STAT), ACTION=TRIM(FACT), ACCESS='STREAM', FORM='UNFORMATTED')
END SUBROUTINE BFOPEN
