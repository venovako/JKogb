MODULE TIMER
  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER FUNCTION GET_SYS_TIME()
    IMPLICIT NONE

    CALL SYSTEM_CLOCK(COUNT=GET_SYS_TIME)
  END FUNCTION GET_SYS_TIME

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER FUNCTION GET_SYS_TRES()
    IMPLICIT NONE

    CALL SYSTEM_CLOCK(COUNT_RATE=GET_SYS_TRES)
  END FUNCTION GET_SYS_TRES

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE TIMER
