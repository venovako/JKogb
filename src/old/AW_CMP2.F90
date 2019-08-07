!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! sorting by type
  ! same type ==> sorting by magnitude
  !   same magnitude ==> sorting by subdiagonals (bands)
  !     outer bands first
  !       within a band, the lower elements first
  FUNCTION AW_CMP2(PA, PB)
    IMPLICIT NONE
    INTEGER(c_intptr_t), INTENT(IN), VALUE :: PA, PB
    INTEGER(c_int) :: AW_CMP2

    TYPE(AW), POINTER :: A, B
    INTEGER :: AB, BB, AP, BP, AQ, BQ

    AW_CMP2 = 0_c_int
    IF (PA .EQ. PB) RETURN

    CALL C_F_POINTER(TRANSFER(PA, C_NULL_PTR), A)
    CALL C_F_POINTER(TRANSFER(PB, C_NULL_PTR), B)

    IF (A%B .LT. 0) THEN
       IF (B%B .GE. 0) THEN
          AW_CMP2 = 1_c_int
          RETURN
       ELSE ! B%B < 0
          BB = -B%B
       END IF
       AB = -A%B
    ELSE ! A%B >= 0
       IF (B%B .LT. 0) THEN
          AW_CMP2 = -1_c_int
          RETURN
       ELSE ! B%B >= 0
          BB = B%B
       END IF
       AB = A%B
    END IF

    IF (A%W .LT. B%W) THEN
       AW_CMP2 = 2_c_int
    ELSE IF (A%W .GT. B%W) THEN
       AW_CMP2 = -2_c_int
    ELSE IF (A%W .EQ. B%W) THEN
       IF (AB .LT. BB) THEN
          AW_CMP2 = 3_c_int
       ELSE IF (AB .GT. BB) THEN
          AW_CMP2 = -3_c_int
       ELSE ! AB = BB
          AP = ABS(A%P)
          BP = ABS(B%P)
          IF (AP .LT. BP) THEN
             AW_CMP2 = 4_c_int
          ELSE IF (AP .GT. BP) THEN
             AW_CMP2 = -4_c_int
          ELSE ! AP = BP
             AQ = ABS(A%Q)
             BQ = ABS(B%Q)
             IF (AQ .LT. BQ) THEN
                AW_CMP2 = 5_c_int
             ELSE IF (AQ .GT. BQ) THEN
                AW_CMP2 = -5_c_int
             ELSE ! AQ = BQ
                AW_CMP2 = 0_c_int
             END IF
          END IF
       END IF
    ELSE ! NaN magnitude(s)
       IF (B%W .EQ. B%W) THEN
          AW_CMP2 = 6_c_int
       ELSE ! NaN(B%W)
          AW_CMP2 = -6_c_int
       END IF
    END IF
  END FUNCTION AW_CMP2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
