MODULE utils
  USE params
  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef NDEBUG
  INTEGER FUNCTION OPEN_LOG(FN, S)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: FN
    INTEGER, INTENT(IN) :: S

    CHARACTER(LEN=LEN_TRIM(FN)+12) :: F
    INTEGER :: U, I

    WRITE (F,'(A,A,I11.11)') TRIM(FN), '.', S
    OPEN(NEWUNIT=U, IOSTAT=I, FILE=F, STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
    IF (I .EQ. 0) THEN
       OPEN_LOG = U
    ELSE ! error
       OPEN_LOG = -1
    END IF
  END FUNCTION OPEN_LOG
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE LOGICAL FUNCTION VERIFY_MIN(STRONG)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: STRONG

    IF (STRONG) THEN
       VERIFY_MIN = (&
            (MIN(QUIET_NAN(-1), 0.0_DWP) .EQ. 0.0_DWP) .AND. &
            (MIN(0.0_DWP, QUIET_NAN(-1)) .EQ. 0.0_DWP))
    ELSE ! weak
       VERIFY_MIN = (MIN(QUIET_NAN(-1), 0.0_DWP) .EQ. 0.0_DWP)
    END IF
  END FUNCTION VERIFY_MIN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE LOGICAL FUNCTION VERIFY_MAX(STRONG)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: STRONG

    IF (STRONG) THEN
       VERIFY_MAX = (&
            (MAX(QUIET_NAN(0), -1.0_DWP) .EQ. -1.0_DWP) .AND. &
            (MAX(-1.0_DWP, QUIET_NAN(0)) .EQ. -1.0_DWP))
    ELSE ! weak
       VERIFY_MAX = (MAX(QUIET_NAN(0), -1.0_DWP) .EQ. -1.0_DWP)
    END IF
  END FUNCTION VERIFY_MAX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE LOGICAL FUNCTION VERIFY_MIN_MAX(STRONG)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: STRONG

    VERIFY_MIN_MAX = (VERIFY_MIN(STRONG) .AND. VERIFY_MAX(STRONG))
  END FUNCTION VERIFY_MIN_MAX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE utils
