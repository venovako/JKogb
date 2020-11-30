MODULE timer
  USE params
  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER(KIND=DWP) FUNCTION GET_SYS_TIME()
    IMPLICIT NONE

    CALL SYSTEM_CLOCK(COUNT=GET_SYS_TIME)
  END FUNCTION GET_SYS_TIME

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER(KIND=DWP) FUNCTION GET_SYS_TRES()
    IMPLICIT NONE

    CALL SYSTEM_CLOCK(COUNT_RATE=GET_SYS_TRES)
  END FUNCTION GET_SYS_TRES

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER(KIND=DWP) FUNCTION GET_SYS_TLAP(T)
    IMPLICIT NONE
    INTEGER(KIND=DWP), INTENT(IN) :: T

    INTEGER(KIND=DWP) :: TNOW, TLAP

    TNOW = GET_SYS_TIME()
    TLAP = TNOW - T
    IF (TLAP .LT. 0_DWP) THEN
       CALL SYSTEM_CLOCK(COUNT_MAX=TNOW)
       TLAP = TLAP + TNOW
    END IF
    GET_SYS_TLAP = MAX(TLAP, 0_DWP)
  END FUNCTION GET_SYS_TLAP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE timer
