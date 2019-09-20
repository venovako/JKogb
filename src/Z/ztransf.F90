MODULE ZTRANSF
  USE PARAMS
  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE ZHSVD2D(H, A, U, Z, INFO)
    ! A diagonal
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: H
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    COMPLEX(KIND=DWP) :: V
    REAL(KIND=DWP) :: W

    IF (AIMAG(A(1,1)) .EQ. D_ZERO) THEN
       ! A(1,1) real
       IF (REAL(A(1,1)) .LT. D_ZERO) THEN
          ! A(1,1) negative
          U(1,1) = -U(1,1)
          U(1,2) = -U(1,2)
          A(1,1) = -A(1,1)
       END IF
    ELSE IF (REAL(A(1,1)) .EQ. D_ZERO) THEN
       ! A(1,1) imaginary .NE. 0
       IF (AIMAG(A(1,1)) .LT. D_ZERO) THEN
          ! A(1,1) = i * negative
          ! V = CMPLX(D_ZERO, D_ONE, DWP)
          U(1,1) = CMPLX(-AIMAG(U(1,1)), REAL(U(1,1)), DWP)
          U(1,2) = CMPLX(-AIMAG(U(1,2)), REAL(U(1,2)), DWP)
          A(1,1) = CMPLX(-AIMAG(A(1,1)), REAL(A(1,1)), DWP)
       ELSE
          ! A(1,1) = i * positive
          ! V = CMPLX(D_ZERO, D_MONE, DWP)
          U(1,1) = CMPLX(AIMAG(U(1,1)), -REAL(U(1,1)), DWP)
          U(1,2) = CMPLX(AIMAG(U(1,2)), -REAL(U(1,2)), DWP)
          A(1,1) = CMPLX(AIMAG(A(1,1)), -REAL(A(1,1)), DWP)
       END IF
    ELSE
       ! A(1,1) complex .NE. 0
       W = ABS(A(1,1))
       V = CONJG(A(1,1) / W)
       U(1,1) = V * U(1,1)
       U(1,2) = V * U(1,2)
       A(1,1) = W
    END IF

    IF (AIMAG(A(2,2)) .EQ. D_ZERO) THEN
       ! A(2,2) real
       IF (REAL(A(2,2)) .LT. D_ZERO) THEN
          ! A(2,2) negative
          U(2,1) = -U(2,1)
          U(2,2) = -U(2,2)
          A(2,2) = -A(2,2)
       END IF
    ELSE IF (REAL(A(2,2)) .EQ. D_ZERO) THEN
       ! A(2,2) imaginary .NE. 0
       IF (AIMAG(A(2,2)) .LT. D_ZERO) THEN
          ! A(2,2) = i * negative
          ! V = CMPLX(D_ZERO, D_ONE, DWP)
          U(2,1) = CMPLX(-AIMAG(U(2,1)), REAL(U(2,1)), DWP)
          U(2,2) = CMPLX(-AIMAG(U(2,2)), REAL(U(2,2)), DWP)
          A(2,2) = CMPLX(-AIMAG(A(2,2)), REAL(A(2,2)), DWP)
       ELSE
          ! A(2,2) = i * positive
          ! V = CMPLX(D_ZERO, D_MONE, DWP)
          U(2,1) = CMPLX(AIMAG(U(2,1)), -REAL(U(2,1)), DWP)
          U(2,2) = CMPLX(AIMAG(U(2,2)), -REAL(U(2,2)), DWP)
          A(2,2) = CMPLX(AIMAG(A(2,2)), -REAL(A(2,2)), DWP)
       END IF
    ELSE
       ! A(2,2) complex .NE. 0
       W = ABS(A(2,2))
       V = CONJG(A(2,2) / W)
       U(2,1) = V * U(2,1)
       U(2,2) = V * U(2,2)
       A(2,2) = W
    END IF

    ! ZHSVD2D called
    INFO = 0
  END SUBROUTINE ZHSVD2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE ZHSVD2U(H, A, U, Z, INFO)
    ! A upper triangular, not diagonal
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: H
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    ! ZHSVD2U called
    INFO = 1
  END SUBROUTINE ZHSVD2U

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE ZHSVD2L(H, A, U, Z, INFO)
    ! A lower triangular, not diagonal
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: H
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    ! ZHSVD2L called
    INFO = 2
  END SUBROUTINE ZHSVD2L

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE ZHSVD2G(H, A, U, Z, INFO)
    ! A general
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: H
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    ! ZHSVD2G called
    INFO = 3
  END SUBROUTINE ZHSVD2G

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE ZHSVD2S(H, A, U, Z, INFO)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: H
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    COMPLEX(KIND=DWP) :: V

    ! assume A real, non-negative, diagonal
    ! permute A, U, Z in trigonometric case
    IF ((.NOT. H) .AND. (REAL(A(1,1)) .LT. REAL(A(2,2)))) THEN
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

    IF ((U(1,1) .NE. Z_ONE) .OR. (U(2,1) .NE. Z_ZERO) .OR. (U(1,2) .NE. Z_ZERO) .OR. (U(2,2) .NE. Z_ONE)) INFO = INFO + 4
    IF ((Z(1,1) .NE. Z_ONE) .OR. (Z(2,1) .NE. Z_ZERO) .OR. (Z(1,2) .NE. Z_ZERO) .OR. (Z(2,2) .NE. Z_ONE)) INFO = INFO + 8
  END SUBROUTINE ZHSVD2S

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE ZHSVD2(H, A, U, Z, INFO)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: H
    COMPLEX(KIND=DWP), INTENT(INOUT) :: A(2,2)
    COMPLEX(KIND=DWP), INTENT(OUT) :: U(2,2), Z(2,2)
    INTEGER, INTENT(OUT) :: INFO

    REAL(KIND=DWP) :: W(2,2)

    W(1,1) = ABS(A(1,1))
    W(2,1) = ABS(A(2,1))
    W(1,2) = ABS(A(1,2))
    W(2,2) = ABS(A(2,2))

    IF (.NOT. (W(1,1) .LE. HUGE(D_ZERO))) THEN
       INFO = -1
    ELSE IF (.NOT. (W(2,1) .LE. HUGE(D_ZERO))) THEN
       INFO = -2
    ELSE IF (.NOT. (W(1,2) .LE. HUGE(D_ZERO))) THEN
       INFO = -3
    ELSE IF (.NOT. (W(2,2) .LE. HUGE(D_ZERO))) THEN
       INFO = -4
    ELSE ! A has no NaNs or infinities
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    ! U = I
    U(1,1) = Z_ONE
    U(2,1) = Z_ZERO
    U(1,2) = Z_ZERO
    U(2,2) = Z_ONE

    ! Z = I
    Z(1,1) = Z_ONE
    Z(2,1) = Z_ZERO
    Z(1,2) = Z_ZERO
    Z(2,2) = Z_ONE

    IF (W(2,1) .EQ. D_ZERO) THEN
       ! A upper triangular
       IF (W(1,2) .EQ. D_ZERO) THEN
          ! A diagonal
          CALL ZHSVD2D(H, A, U, Z, INFO)
       ELSE
          ! A upper triangular, not diagonal
          CALL ZHSVD2U(H, A, U, Z, INFO)
       END IF
    ELSE IF (W(1,2) .EQ. D_ZERO) THEN
       ! A lower triangular, not diagonal
       CALL ZHSVD2L(H, A, U, Z, INFO)
    ELSE
       ! A general
       CALL ZHSVD2G(H, A, U, Z, INFO)
    END IF

    CALL ZHSVD2S(H, A, U, Z, INFO)
  END SUBROUTINE ZHSVD2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE ZTRANSF
