MODULE UTILS
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE

  TYPE(c_funptr) :: OldCtrlC = C_NULL_FUNPTR
  INTEGER(c_int), PARAMETER :: SigCtrlC = 2_c_int ! SIGINT
  INTEGER(c_int), VOLATILE :: CtrlC = 0_c_int

  ! INTERFACE
  !    RECURSIVE FUNCTION PAR_SORT(T, A, N, C, X) BIND(C,NAME='par_sort')
  !      USE, INTRINSIC :: ISO_C_BINDING
  !      IMPLICIT NONE
  !      INTEGER(KIND=c_int), INTENT(IN), VALUE :: T
  !      TYPE(c_ptr), INTENT(IN), VALUE :: A, X
  !      INTEGER(KIND=c_size_t), INTENT(IN), VALUE :: N
  !      TYPE(c_funptr), INTENT(IN), VALUE :: C
  !      INTEGER(KIND=c_int) :: PAR_SORT
  !    END FUNCTION PAR_SORT
  ! END INTERFACE

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

  SUBROUTINE OnCtrlC(S) BIND(C)
    IMPLICIT NONE
    INTEGER(c_int), INTENT(IN), VALUE :: S
    IF (S .EQ. SigCtrlC) CtrlC = 1_c_int
  END SUBROUTINE OnCtrlC

  SUBROUTINE SetCtrlC
    IMPLICIT NONE
    INTERFACE
       FUNCTION C_SIGNAL(S, H) BIND(C,NAME='signal')
         USE, INTRINSIC :: ISO_C_BINDING
         INTEGER(c_int), INTENT(IN), VALUE :: S
         TYPE(c_funptr), INTENT(IN), VALUE :: H
         TYPE(c_funptr) :: C_SIGNAL
       END FUNCTION C_SIGNAL
    END INTERFACE

    OldCtrlC = C_SIGNAL(SigCtrlC, C_FUNLOC(OnCtrlC))
  END SUBROUTINE SetCtrlC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE UTILS
