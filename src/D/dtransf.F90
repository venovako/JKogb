MODULE DTRANSF
  USE PARAMS
  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE DHSVD2D(H, A, U, Z, INFO)
    ! A diagonal
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: H
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    REAL(KIND=DWP) :: W

    IF (A(1,1) .LT. D_ZERO) THEN
       ! A(1,1) negative
       U(1,1) = -U(1,1)
       U(1,2) = -U(1,2)
       A(1,1) = -A(1,1)
    END IF

    IF (A(2,2) .LT. D_ZERO) THEN
       ! A(2,2) negative
       U(2,1) = -U(2,1)
       U(2,2) = -U(2,2)
       A(2,2) = -A(2,2)
    END IF

    ! DHSVD2D called
    INFO = 0
  END SUBROUTINE DHSVD2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE DHSVD2U(H, A, U, Z, INFO)
    ! A upper triangular, not diagonal
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: H
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    ! DHSVD2U called
    INFO = 1
  END SUBROUTINE DHSVD2U

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE DHSVD2L(H, A, U, Z, INFO)
    ! A lower triangular, not diagonal
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: H
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    ! DHSVD2L called
    INFO = 2
  END SUBROUTINE DHSVD2L

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE DHSVD2G(H, A, U, Z, INFO)
    ! A general
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: H
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    ! DHSVD2G called
    INFO = 3
  END SUBROUTINE DHSVD2G

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE DHSVD2S(H, A, U, Z, INFO)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: H
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    REAL(KIND=DWP) :: V

    ! assume A real, non-negative, diagonal
    ! permute A, U, Z in trigonometric case
    IF ((.NOT. H) .AND. (A(1,1) .LT. A(2,2))) THEN
       ! swap the rows of U
       V = U(1,1)
       U(1,1) = U(2,1)
       U(2,1) = V
       V = U(1,2)
       U(1,2) = U(2,2)
       U(2,2) = V
       ! swap the diagonal elements of A
       V = A(1,1)
       A(1,1) = A(2,2)
       A(2,2) = V
       ! swap the columns of Z
       V = Z(1,1)
       Z(1,1) = Z(1,2)
       Z(1,2) = V
       V = Z(2,1)
       Z(2,1) = Z(2,2)
       Z(2,2) = V
    END IF

    IF ((U(1,1) .NE. D_ONE) .OR. (U(2,1) .NE. D_ZERO) .OR. (U(1,2) .NE. D_ZERO) .OR. (U(2,2) .NE. D_ONE)) INFO = INFO + 4
    IF ((Z(1,1) .NE. D_ONE) .OR. (Z(2,1) .NE. D_ZERO) .OR. (Z(1,2) .NE. D_ZERO) .OR. (Z(2,2) .NE. D_ONE)) INFO = INFO + 8
  END SUBROUTINE DHSVD2S

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE DHSVD2(H, A, U, Z, INFO)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: H
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2)
    REAL(KIND=DWP), INTENT(OUT) :: U(2,2), Z(2,2)
    INTEGER, INTENT(OUT) :: INFO

    IF (.NOT. (ABS(A(1,1)) .LE. HUGE(D_ZERO))) THEN
       INFO = -1
    ELSE IF (.NOT. (ABS(A(2,1)) .LE. HUGE(D_ZERO))) THEN
       INFO = -2
    ELSE IF (.NOT. (ABS(A(1,2)) .LE. HUGE(D_ZERO))) THEN
       INFO = -3
    ELSE IF (.NOT. (ABS(A(2,2)) .LE. HUGE(D_ZERO))) THEN
       INFO = -4
    ELSE ! A has no NaNs or infinities
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    ! U = I
    U(1,1) = D_ONE
    U(2,1) = D_ZERO
    U(1,2) = D_ZERO
    U(2,2) = D_ONE

    ! Z = I
    Z(1,1) = D_ONE
    Z(2,1) = D_ZERO
    Z(1,2) = D_ZERO
    Z(2,2) = D_ONE

    IF (A(2,1) .EQ. D_ZERO) THEN
       ! A upper triangular
       IF (A(1,2) .EQ. D_ZERO) THEN
          ! A diagonal
          CALL DHSVD2D(H, A, U, Z, INFO)
       ELSE
          ! A upper triangular, not diagonal
          CALL DHSVD2U(H, A, U, Z, INFO)
       END IF
    ELSE IF (A(1,2) .EQ. D_ZERO) THEN
       ! A lower triangular, not diagonal
       CALL DHSVD2L(H, A, U, Z, INFO)
    ELSE
       ! A general
       CALL DHSVD2G(H, A, U, Z, INFO)
    END IF

    CALL DHSVD2S(H, A, U, Z, INFO)
  END SUBROUTINE DHSVD2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE DTRANSF
