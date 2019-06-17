MODULE UTILS
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE

  INTERFACE
     SUBROUTINE PAR_SORT(A, N, C) BIND(C,NAME='par_sort')
       USE, INTRINSIC :: ISO_C_BINDING
       TYPE(c_ptr), INTENT(IN), VALUE :: A
       INTEGER(c_size_t), INTENT(IN), VALUE :: N
       TYPE(c_funptr), INTENT(IN), VALUE :: C
     END SUBROUTINE PAR_SORT
  END INTERFACE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER FUNCTION OPEN_LOG(FN)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: FN

    INTEGER, SAVE :: S = 0
    CHARACTER(LEN=LEN_TRIM(FN)+12) :: F
    INTEGER :: U, I

    WRITE (F,'(A,A,I11.11)') TRIM(FN), '.', S
    OPEN(NEWUNIT=U, IOSTAT=I, FILE=F, STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
    IF (I .EQ. 0) THEN
       OPEN_LOG = U
       S = S + 1
    ELSE ! error
       OPEN_LOG = -1
    END IF
  END FUNCTION OPEN_LOG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE UTILS
