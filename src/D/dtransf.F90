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

    ! INFO will stay 0 iff no transformations have been applied
    INFO = 0

    IF (SIGN(D_ONE, A(1,1)) .EQ. D_MONE) THEN
       ! A(1,1) negative
       U(1,1) = -U(1,1)
       U(1,2) = -U(1,2)
       A(1,1) = -A(1,1)
       INFO = 1
    END IF

    IF (SIGN(D_ONE, A(2,2)) .EQ. D_MONE) THEN
       ! A(2,2) negative
       U(2,1) = -U(2,1)
       U(2,2) = -U(2,2)
       A(2,2) = -A(2,2)
       INFO = INFO + 2
    END IF

    ! DHSVD2D called
    IF (INFO .NE. 0) INFO = 1
  END SUBROUTINE DHSVD2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DHSVD2U(H, A, U, Z, INFO)
    ! A upper triangular, not diagonal
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: H
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    REAL(KIND=DWP) :: V(2,2), W(2,2), B(2,2), CS, SN, TN, CH, SH, TH

    EXTERNAL :: DGEMM

    INFO = 0

    ! TODO: transform
    IF (H) THEN
       IF (A(2,1) .NE. D_ZERO) THEN
          CONTINUE
       ELSE IF (ABS(A(1,2)) .LT. ABS(A(1,1))) THEN
          CONTINUE
       ELSE ! (A(2,1) .EQ. 0) .AND. (ABS(A(1,1)) .EQ. ABS(A(1,2)))
          ! |tanh| .GE. 1
          INFO=-6
          RETURN
       END IF
    ELSE
       IF (A(2,1) .NE. D_ZERO) THEN
          CONTINUE
       ELSE ! A(2,1) .EQ. 0
          CONTINUE
       END IF
    END IF

    CALL DGEMM('N', 'N', 2, 2, 2, D_ONE, V, 2, U, 2, D_ZERO, B, 2)
    U = B
    CALL DGEMM('N', 'N', 2, 2, 2, D_ONE, V, 2, A, 2, D_ZERO, B, 2)
    A = B
    CALL DGEMM('N', 'N', 2, 2, 2, D_ONE, Z, 2, W, 2, D_ZERO, B, 2)
    Z = B
    CALL DGEMM('N', 'N', 2, 2, 2, D_ONE, A, 2, W, 2, D_ZERO, B, 2)
    A = B

    A(2,1) = D_ZERO
    A(1,2) = D_ZERO
    CALL DHSVD2D(H, A, U, Z, INFO)

    ! DHSVD2U called
    INFO = 2
  END SUBROUTINE DHSVD2U

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DHSVD2T(H, A, U, Z, INFO)
    ! A upper antitriangular, not antidiagonal
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: H
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    REAL(KIND=DWP) :: V(2,2), W(2,2), B(2,2), CS, SN, TN, CH, SH, TH

    EXTERNAL :: DGEMM

    IF (.NOT. H) THEN
       INFO = -HUGE(0)
       RETURN
    END IF
    INFO = 0

    ! TODO: transform
    IF (A(2,1) .NE. D_ZERO) THEN
       CONTINUE
    ELSE IF (ABS(A(1,1)) .LT. ABS(A(1,2))) THEN
       V(1,1) = D_ZERO
       V(2,1) = D_ONE
       V(1,2) = D_MONE
       V(2,2) = D_ZERO
       TH = -A(1,1) / A(1,2)
       !DIR$ FMA
       CH = D_ONE / SQRT(D_ONE - TH * TH)
       SH = TH * CH
       W(1,1) = CH
       W(2,1) = SH
       W(1,2) = SH
       W(2,2) = CH
    ELSE ! (A(2,1) .EQ. 0) .AND. (ABS(A(1,1)) .EQ. ABS(A(1,2)))
       ! |tanh| .GE. 1
       INFO = -7
       RETURN
    END IF

    CALL DGEMM('N', 'N', 2, 2, 2, D_ONE, V, 2, U, 2, D_ZERO, B, 2)
    U = B
    CALL DGEMM('N', 'N', 2, 2, 2, D_ONE, V, 2, A, 2, D_ZERO, B, 2)
    A = B
    CALL DGEMM('N', 'N', 2, 2, 2, D_ONE, Z, 2, W, 2, D_ZERO, B, 2)
    Z = B
    CALL DGEMM('N', 'N', 2, 2, 2, D_ONE, A, 2, W, 2, D_ZERO, B, 2)
    A = B

    A(2,1) = D_ZERO
    A(1,2) = D_ZERO
    CALL DHSVD2D(H, A, U, Z, INFO)

    ! DHSVD2T called
    INFO = 3
  END SUBROUTINE DHSVD2T

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DHSVD2G(H, A, U, Z, INFO)
    ! A general
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: H
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    REAL(KIND=DWP) :: C, S, R

    REAL(KIND=DWP), EXTERNAL :: DNRM2
    EXTERNAL :: DLARTG, DROT, DSWAP

    IF (DNRM2(2, A(1,1), 1) .LT. DNRM2(2, A(1,2), 1)) THEN
       ! column swap of A
       CALL DSWAP(2, A(1,1), 1, A(1,2), 1)
       ! column swap of Z
       IF (.NOT. H) CALL DSWAP(2, Z(1,1), 1, Z(1,2), 1)
       ! record the swap
       INFO = 1
    ELSE
       INFO = 0
    END IF

    IF (ABS(A(1,1)) .LT. ABS(A(2,1))) THEN
       ! row swap of U
       CALL DSWAP(2, U(1,1), 2, U(2,1), 2)
       ! row swap of A
       CALL DSWAP(2, A(1,1), 2, A(2,1), 2)
    END IF

    ! QR factorization of A
    CALL DLARTG(A(1,1), A(2,1), C, S, R)
    IF (.NOT. (ABS(R) .LE. HUGE(D_ZERO))) THEN
       INFO = -5
       RETURN
    END IF
    A(1,1) = R
    A(2,1) = D_ZERO
    CALL DROT(1, A(1,2), 2, A(2,2), 2, C, S)
    ! premultiply U by Q^T
    CALL DROT(2, U(1,1), 2, U(2,1), 2, C, S)

    IF (A(1,2) .EQ. D_ZERO) THEN
       ! A diagonal
       IF (H .AND. (INFO .EQ. 1)) THEN
          ! swap the rows of U
          CALL DSWAP(2, U(1,1), 2, U(2,1), 2)
          ! swap the diagonal elements of A
          CALL DSWAP(1, A(1,1), 1, A(2,2), 1)
       END IF
       CALL DHSVD2D(H, A, U, Z, INFO)
    ELSE IF (H .AND. (INFO .EQ. 1)) THEN
       ! column swap of A
       CALL DSWAP(2, A(1,1), 1, A(1,2), 1)
       ! A upper antitriangular, not antidiagonal (X .NE. 0)
       !     | X R | <- R .NE. 0 the largest
       ! A = | x 0 |    element by magnitude
       CALL DHSVD2T(H, A, U, Z, INFO)
    ELSE
       ! A upper triangular, not diagonal (X .NE. 0)
       !     | R X | <- R .NE. 0 the largest
       ! A = | 0 x |    element by magnitude
       CALL DHSVD2U(H, A, U, Z, INFO)
    END IF
  END SUBROUTINE DHSVD2G

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DHSVD2S(H, A, U, Z, INFO)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: H
    REAL(KIND=DWP), INTENT(INOUT) :: A(2,2), U(2,2), Z(2,2)
    INTEGER, INTENT(INOUT) :: INFO

    EXTERNAL :: DSWAP

    ! assume A real, non-negative, diagonal
    ! permute A, U, Z in trigonometric case
    IF ((.NOT. H) .AND. (A(1,1) .LT. A(2,2))) THEN
       ! swap the rows of U
       CALL DSWAP(2, U(1,1), 2, U(2,1), 2)
       ! swap the diagonal elements of A
       CALL DSWAP(1, A(1,1), 1, A(2,2), 1)
       ! swap the columns of Z
       CALL DSWAP(2, Z(1,1), 1, Z(1,2), 1)
    END IF

    ! check if U is identity and record in INFO if it is not
    IF ((U(1,1) .NE. D_ONE) .OR. (U(2,1) .NE. D_ZERO) .OR. (U(1,2) .NE. D_ZERO) .OR. (U(2,2) .NE. D_ONE)) INFO = INFO + 4
    ! check if Z is identity and record in INFO if it is not
    IF ((Z(1,1) .NE. D_ONE) .OR. (Z(2,1) .NE. D_ZERO) .OR. (Z(1,2) .NE. D_ZERO) .OR. (Z(2,2) .NE. D_ONE)) INFO = INFO + 8
  END SUBROUTINE DHSVD2S

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DHSVD2(H, A, U, Z, INFO)
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

    IF ((A(2,1) .EQ. D_ZERO) .AND. (A(1,2) .EQ. D_ZERO)) THEN
       ! A diagonal
       CALL DHSVD2D(H, A, U, Z, INFO)
    ELSE
       ! A general
       CALL DHSVD2G(H, A, U, Z, INFO)
    END IF
    IF (INFO .LT. 0) RETURN

    CALL DHSVD2S(H, A, U, Z, INFO)
  END SUBROUTINE DHSVD2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE DTRANSF
