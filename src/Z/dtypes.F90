MODULE DTYPES
  USE PARAMS
  USE VN_SORT_F
  IMPLICIT NONE

  ABSTRACT INTERFACE
     PURE FUNCTION AMAG(APP, AQP, APQ, AQQ, JP, JQ)
       USE PARAMS
       IMPLICIT NONE
       COMPLEX(KIND=DWP), INTENT(IN) :: APP, AQP, APQ, AQQ
       INTEGER, INTENT(IN) :: JP, JQ
       REAL(KIND=DWP) :: AMAG
     END FUNCTION AMAG
  END INTERFACE

  TYPE AMP
     PROCEDURE(AMAG), POINTER, NOPASS :: AM
  END TYPE AMP

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

  SUBROUTINE DZBW_OUT(OU, DZ)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: OU
    TYPE(DZBW), INTENT(IN) :: DZ

    IF (OU .GE. 0) THEN
       WRITE (UNIT=OU,FMT=1) '{ ', DZ%W, ', ', DZ%P, ', ', DZ%Q, ', ', DZ%B, ' }'
    ELSE ! OU < 0
       WRITE (UNIT=-OU,FMT=1,ADVANCE='NO') '{ ', DZ%W, ', ', DZ%P, ', ', DZ%Q, ', ', DZ%B, ' }'
    END IF
1   FORMAT(A,ES25.17E3,A,I11,A,I11,A,I11,A)
  END SUBROUTINE DZBW_OUT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE SUBROUTINE DZBW_GEN(APP, AQP, APQ, AQQ, JP, JQ, P, Q, AM, DZ)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(IN) :: APP, AQP, APQ, AQQ
    INTEGER, INTENT(IN) :: JP, JQ, P, Q
    PROCEDURE(AMAG) :: AM
    TYPE(DZBW), INTENT(OUT) :: DZ

    DZ%W = AM(APP, AQP, APQ, AQQ, JP, JQ)
    IF (JP .GE. 0) THEN
       IF (JQ .GE. 0) THEN
          DZ%P = P
          DZ%Q = Q
          DZ%B = Q - P
       ELSE ! JQ < 0
          DZ%P = P
          DZ%Q = -Q
          DZ%B = P - Q
       END IF
    ELSE ! JP < 0
       IF (JQ .GE. 0) THEN
          DZ%P = -P
          DZ%Q = Q
          DZ%B = P - Q
       ELSE ! JQ < 0
          DZ%P = -P
          DZ%Q = -Q
          DZ%B = Q - P
       END IF
    END IF
  END SUBROUTINE DZBW_GEN

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

    DZBW_CMP = 0_c_int
    IF (PA .EQ. PB) RETURN

    CALL C_F_POINTER(TRANSFER(PA, C_NULL_PTR), A)
    CALL C_F_POINTER(TRANSFER(PB, C_NULL_PTR), B)

    IF (A%B .LT. 0) THEN
       IF (B%B .GE. 0) THEN
          DZBW_CMP = 1_c_int
          RETURN
       ELSE ! B%B < 0
          BB = -B%B
       END IF
       AB = -A%B
    ELSE ! A%B >= 0
       IF (B%B .LT. 0) THEN
          DZBW_CMP = -1_c_int
          RETURN
       ELSE ! B%B >= 0
          BB = B%B
       END IF
       AB = A%B
    END IF

    IF (A%W .LT. B%W) THEN
       DZBW_CMP = 2_c_int
    ELSE IF (A%W .GT. B%W) THEN
       DZBW_CMP = -2_c_int
    ELSE IF (A%W .EQ. B%W) THEN
       IF (AB .LT. BB) THEN
          DZBW_CMP = 3_c_int
       ELSE IF (AB .GT. BB) THEN
          DZBW_CMP = -3_c_int
       ELSE ! AB = BB
          AP = ABS(A%P)
          BP = ABS(B%P)
          IF (AP .LT. BP) THEN
             DZBW_CMP = 4_c_int
          ELSE IF (AP .GT. BP) THEN
             DZBW_CMP = -4_c_int
          ELSE ! AP = BP
             AQ = ABS(A%Q)
             BQ = ABS(B%Q)
             IF (AQ .LT. BQ) THEN
                DZBW_CMP = 5_c_int
             ELSE IF (AQ .GT. BQ) THEN
                DZBW_CMP = -5_c_int
             ELSE ! AQ = BQ
                DZBW_CMP = 0_c_int
             END IF
          END IF
       END IF
    ELSE ! NaN magnitude(s)
       IF (B%W .EQ. B%W) THEN
          DZBW_CMP = 6_c_int
       ELSE ! NaN(B%W)
          DZBW_CMP = -6_c_int
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

  PURE SUBROUTINE DZBW_NCP(NN, A, N_2, S, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NN, N_2
    TYPE(DZBW), INTENT(IN) :: A(NN)
    INTEGER, INTENT(OUT) :: S(N_2), INFO

    INTEGER :: I, J, K, AP, AQ, BP, BQ
    LOGICAL :: C

    IF (NN .LT. 0) THEN
       INFO = -1
    ELSE IF (N_2 .LT. 0) THEN
       INFO = -3
    ELSE IF (N_2 .GT. NN) THEN
       INFO = -3
    ELSE
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN
    IF (NN .EQ. 0) RETURN
    IF (N_2 .EQ. 0) RETURN

    J = 1
    DO I = 1, N_2
       S(I) = J
       J = J + 1
       DO WHILE (J .LE. NN)
          AP = ABS(A(J)%P)
          AQ = ABS(A(J)%Q)
          C = .FALSE.
          DO K = I, 1, -1
             BP = ABS(A(S(K))%P)
             BQ = ABS(A(S(K))%Q)
             IF ((AP .EQ. BP) .OR. (AP .EQ. BQ) .OR. (AQ .EQ. BP) .OR. (AQ .EQ. BQ)) THEN
                C = .TRUE.
                EXIT
             END IF
          END DO
          IF (C) THEN
             J = J + 1
          ELSE
             EXIT
          END IF
       END DO
       INFO = I
       IF (J .GT. NN) EXIT
    END DO

    !DIR$ VECTOR ALWAYS
    DO I = INFO+1, N_2
       S(I) = 0
    END DO
  END SUBROUTINE DZBW_NCP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE DTYPES
