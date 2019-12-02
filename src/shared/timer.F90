MODULE TIMER
  USE PARAMS
  USE VN_TIMER_F
  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER FUNCTION GET_THREAD_NS()
    IMPLICIT NONE
    GET_THREAD_NS = INT(VN_GET_THREAD_NS())
  END FUNCTION GET_THREAD_NS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER FUNCTION GET_SYS_US()
    IMPLICIT NONE
    GET_SYS_US = INT(VN_GET_SYS_US())
  END FUNCTION GET_SYS_US

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  REAL(KIND=c_double) FUNCTION TIMER2DBLE(CLK)
    IMPLICIT NONE
    INTEGER(KIND=c_int64_t), INTENT(IN), TARGET :: CLK(4)
#ifdef USE_TSC
    INTEGER(KIND=c_int64_t), TARGET :: SEC, REM
    TIMER2DBLE = TSC_LAP(CLK(3), CLK(1), CLK(2), SEC, REM)
#else
    TIMER2DBLE = (CLK(2) - CLK(1)) / REAL(CLK(3), c_double)
#endif
  END FUNCTION TIMER2DBLE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE TIMER_INIT(CLK)
    IMPLICIT NONE
    INTEGER(KIND=c_int64_t), INTENT(OUT), TARGET :: CLK(4)
#ifdef USE_TSC
    CLK(3) = TSC_GET_FREQ_HZ(CLK(4))
#else
    CALL SYSTEM_CLOCK(CLK(2), CLK(3), CLK(4))
#endif
    CLK(1) = 0_c_int64_t
    CLK(2) = 0_c_int64_t
  END SUBROUTINE TIMER_INIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE TIMER_START(CLK)
    IMPLICIT NONE
    INTEGER(KIND=c_int64_t), INTENT(OUT), TARGET :: CLK(4)
#ifdef USE_TSC
    INTEGER(KIND=c_int32_t), TARGET :: AUX
    CLK(1) = INT(RDTSC_BEG(AUX))
#else
    CALL SYSTEM_CLOCK(CLK(1))
#endif
  END SUBROUTINE TIMER_START

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE TIMER_STOP(CLK)
    IMPLICIT NONE
    INTEGER(KIND=c_int64_t), INTENT(OUT), TARGET :: CLK(4)
#ifdef USE_TSC
    INTEGER(KIND=c_int32_t), TARGET :: AUX
    CLK(2) = INT(RDTSC_END(AUX))
#else
    CALL SYSTEM_CLOCK(CLK(2))
#endif
  END SUBROUTINE TIMER_STOP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE TIMER
