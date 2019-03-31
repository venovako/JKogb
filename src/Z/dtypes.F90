MODULE DTYPES
  USE PARAMS
  USE VN_SORT_F
  IMPLICIT NONE

  TYPE, BIND(C) :: DZBW
     ! W = |A_pq|+|A_qp|
     REAL(KIND=DWP) :: W ! weight
     INTEGER :: P ! row * sign(J_p)
     INTEGER :: Q ! column * sign(J_q)
     ! |B| = |Q| - |P| > 0
     INTEGER :: B ! band * sign(J_p) * sign(J_q)
  END TYPE DZBW

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! sorting by type
  ! same type ==> sorting by magnitude
  !   same magnitude ==> sorting by subdiagonals (bands)
  !     outer bands first
  !       within a band, the lower elements first
  FUNCTION DZBW_CMP(PA, PB)
    IMPLICIT NONE
    INTEGER(c_intptr_t), INTENT(IN), VALUE :: PA, PB
    INTEGER(c_int) :: DZBW_CMP

    TYPE(DZBW), POINTER :: A, B
    INTEGER :: AB, BB, AP, BP, AQ, BQ

    DZBW_CMP = 0
    IF (PA .EQ. PB) RETURN

    CALL C_F_POINTER(TRANSFER(PA, C_NULL_PTR), A)
    CALL C_F_POINTER(TRANSFER(PB, C_NULL_PTR), B)

    IF (A%B .LT. 0) THEN
       IF (B%B .GE. 0) THEN
          DZBW_CMP = 1
          RETURN
       ELSE ! B%B < 0
          BB = -B%B
       END IF
       AB = -A%B
    ELSE ! A%B >= 0
       IF (B%B .LT. 0) THEN
          DZBW_CMP = -1
          RETURN
       ELSE ! B%B >= 0
          BB = B%B
       END IF
       AB = A%B
    END IF

    IF (A%W .LT. B%W) THEN
       DZBW_CMP = 2
    ELSE IF (A%W .GT. B%W) THEN
       DZBW_CMP = -2
    ELSE IF (A%W .EQ. B%W) THEN
       IF (AB .LT. BB) THEN
          DZBW_CMP = 3
       ELSE IF (AB .GT. BB) THEN
          DZBW_CMP = -3
       ELSE ! AB = BB
          AP = ABS(A%P)
          BP = ABS(B%P)
          IF (AP .LT. BP) THEN
             DZBW_CMP = 4
          ELSE IF (AP .GT. BP) THEN
             DZBW_CMP = -4
          ELSE ! AP = BP
             AQ = ABS(A%Q)
             BQ = ABS(B%Q)
             IF (AQ .LT. BQ) THEN
                DZBW_CMP = 5
             ELSE IF (AQ .GT. BQ) THEN
                DZBW_CMP = -5
             ELSE ! AQ = BQ
                DZBW_CMP = 0
             END IF
          END IF
       END IF
    ELSE ! NaN magnitude(s)
       IF (B%W .EQ. B%W) THEN
          DZBW_CMP = 6
       ELSE ! NaN(B%W)
          DZBW_CMP = -6
       END IF
    END IF
  END FUNCTION DZBW_CMP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DZBW_SRT(NN, A, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NN
    TYPE(DZBW), INTENT(INOUT), TARGET :: A(NN)
    INTEGER, INTENT(OUT) :: INFO

    IF (NN .LT. 0) THEN
       INFO = -1
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN
    IF (NN .EQ. 0) RETURN

    CALL VN_QSORT(C_LOC(A), INT(NN,c_size_t), C_SIZEOF(A(1)), C_FUNLOC(DZBW_CMP))
  END SUBROUTINE DZBW_SRT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DZBW_NCP(NN, A, S, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NN
    TYPE(DZBW), INTENT(IN) :: A(NN)
    INTEGER, INTENT(OUT) :: S(*), INFO

    IF (NN .LT. 0) THEN
       INFO = -1
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN
    IF (NN .EQ. 0) RETURN
  END SUBROUTINE DZBW_NCP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE DTYPES
