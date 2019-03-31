MODULE DTYPES
  USE PARAMS
  USE VN_SORT_F
  IMPLICIT NONE

  TYPE, BIND(C) :: DZBW
     ! W = |A_pq|+|A_qp|
     REAL(KIND=DWP) :: W ! weight
     INTEGER :: P ! row
     INTEGER :: Q ! column
     ! B = Q - P > 0
     INTEGER :: B ! band
  END TYPE DZBW

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! sorting by magnitude
  ! same magnitude ==> sorting by subdiagonals (bands)
  !   outer bands first
  !     within a band, the lower elements first
  FUNCTION DZBW_CMP(PA, PB)
    IMPLICIT NONE
    INTEGER(c_intptr_t), INTENT(IN), VALUE :: PA, PB
    INTEGER(c_int) :: DZBW_CMP

    TYPE(DZBW), POINTER :: A, B

    DZBW_CMP = 0
    IF (PA .EQ. PB) RETURN

    CALL C_F_POINTER(TRANSFER(PA, C_NULL_PTR), A)
    CALL C_F_POINTER(TRANSFER(PB, C_NULL_PTR), B)

    IF (A%W .LT. B%W) THEN
       DZBW_CMP = 1
    ELSE IF (A%W .GT. B%W) THEN
       DZBW_CMP = -1
    ELSE IF (A%W .EQ. B%W) THEN
       IF (A%B .LT. B%B) THEN
          DZBW_CMP = 2
       ELSE IF (A%B .GT. B%B) THEN
          DZBW_CMP = -2
       ELSE IF (A%P .LT. B%P) THEN
          DZBW_CMP = 3
       ELSE IF (A%P .GT. B%P) THEN
          DZBW_CMP = -3
       ELSE ! all equal
          DZBW_CMP = 0
       END IF
    ELSE ! NaN magnitude(s)
       IF (B%W .EQ. B%W) THEN
          DZBW_CMP = 4
       ELSE ! NaN(B%W)
          DZBW_CMP = -4
       END IF
    END IF
  END FUNCTION DZBW_CMP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DZBW_SRT(N, A, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    TYPE(DZBW), INTENT(INOUT), TARGET :: A(N)
    INTEGER, INTENT(OUT) :: INFO

    IF (N .LT. 0) THEN
       INFO = -1
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN
    IF (N .EQ. 0) RETURN

    CALL VN_QSORT(C_LOC(A), INT(N,c_size_t), C_SIZEOF(A(1)), C_FUNLOC(DZBW_CMP))
  END SUBROUTINE DZBW_SRT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE DTYPES
