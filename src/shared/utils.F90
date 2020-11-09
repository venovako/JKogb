MODULE UTILS
  USE, INTRINSIC :: ISO_C_BINDING
#ifndef USE_NVIDIA
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: ATOMIC_INT_KIND
#endif
  IMPLICIT NONE

#ifdef USE_NVIDIA
  INTEGER, PARAMETER :: ATOMIC_INT_KIND = c_int
#endif
  INTEGER(KIND=ATOMIC_INT_KIND), VOLATILE, PRIVATE :: CtrlC = 0_ATOMIC_INT_KIND
  INTEGER(KIND=c_int), PARAMETER, PRIVATE :: SigCtrlC = 2_c_int ! SIGINT
  TYPE(c_funptr), PRIVATE :: OldCtrlC = C_NULL_FUNPTR

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  LOGICAL FUNCTION IsCtrlC()
    IMPLICIT NONE

    IF (CtrlC .EQ. 0_ATOMIC_INT_KIND) THEN
       IsCtrlC = .FALSE.
    ELSE ! CtrlC
       IsCtrlC = .TRUE.
    END IF
  END FUNCTION IsCtrlC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE OnCtrlC(S) BIND(C)
    IMPLICIT NONE
    INTEGER(KIND=c_int), INTENT(IN), VALUE :: S

    IF (S .EQ. SigCtrlC) CtrlC = 1_ATOMIC_INT_KIND
  END SUBROUTINE OnCtrlC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE SetCtrlC
    IMPLICIT NONE
    INTERFACE
       FUNCTION C_SIGNAL(S, H) BIND(C,NAME='signal')
         USE, INTRINSIC :: ISO_C_BINDING
         INTEGER(KIND=c_int), INTENT(IN), VALUE :: S
         TYPE(c_funptr), INTENT(IN), VALUE :: H
         TYPE(c_funptr) :: C_SIGNAL
       END FUNCTION C_SIGNAL
    END INTERFACE
    OldCtrlC = C_SIGNAL(SigCtrlC, C_FUNLOC(OnCtrlC))
  END SUBROUTINE SetCtrlC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE UTILS
