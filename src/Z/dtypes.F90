MODULE DTYPES
  USE TIMER
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

  ABSTRACT INTERFACE
     PURE FUNCTION ACVG(APP, AQP, APQ, AQQ, JP, JQ)
       USE PARAMS
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: JP, JQ
       COMPLEX(KIND=DWP), INTENT(IN) :: APP, AQP, APQ, AQQ
       INTEGER :: ACVG
     END FUNCTION ACVG
  END INTERFACE

  ABSTRACT INTERFACE
     PURE SUBROUTINE ATRU(N, P, Q, NN, INFO)
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: N, NN
       INTEGER, INTENT(OUT) :: P(NN), Q(NN), INFO
     END SUBROUTINE ATRU
  END INTERFACE

  TYPE APROC
     PROCEDURE(AMAG), POINTER, NOPASS :: MAG
     PROCEDURE(VN_QSORT_CMP), POINTER, NOPASS :: CMP
     PROCEDURE(ACVG), POINTER, NOPASS :: CVG
     PROCEDURE(ATRU), POINTER, NOPASS :: TRU
  END TYPE APROC

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

  SUBROUTINE DZBW_OUT(OU, HDR, NN, DZ, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: OU, NN, SL, STEP(SL)
    CHARACTER(LEN=*), INTENT(IN) :: HDR
    TYPE(DZBW), INTENT(IN) :: DZ(NN)
    INTEGER, INTENT(OUT) :: INFO

    INTEGER :: I, J

    IF (OU .EQ. -1) THEN
       INFO = -1
    ELSE IF (NN .LT. 0) THEN
       INFO = -3
    ELSE IF (SL .LT. 0) THEN
       INFO = -5
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    INFO = GET_THREAD_NS()

    IF (LEN_TRIM(HDR) .GT. 0) THEN
       WRITE (OU,'(A)') TRIM(HDR)
    ELSE IF (SL .EQ. 0) THEN
       WRITE (OU,'(A)') '"J","W","P","Q","B"'
    ELSE ! SL > 0
       WRITE (OU,'(A)') '"I","J","W","P","Q","B"'
    END IF

    IF (SL .EQ. 0) THEN
       DO J = 1, NN
          WRITE (OU,'(I11,A,ES25.17E3,3(A,I11))') J, ',', DZ(J)%W, ',', DZ(J)%P, ',', DZ(J)%Q, ',', DZ(J)%B
       END DO
    ELSE ! SL > 0
       DO I = 1, SL
          J = STEP(I)
          WRITE (OU,'(2(I11,A),ES25.17E3,3(A,I11))') I, ',', J, ',', DZ(J)%W, ',', DZ(J)%P, ',', DZ(J)%Q, ',', DZ(J)%B
       END DO
    END IF

    INFO = GET_THREAD_NS() - INFO
  END SUBROUTINE DZBW_OUT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! sorting by magnitude
  !   same magnitude ==> sorting by subdiagonals (bands)
  !     outer bands first
  !       within a band, the lower elements first
  FUNCTION DZBW_CMP1(PA, PB)
    IMPLICIT NONE
    INTEGER(c_intptr_t), INTENT(IN), VALUE :: PA, PB
    INTEGER(c_int) :: DZBW_CMP1

    TYPE(DZBW), POINTER :: A, B
    INTEGER :: AB, BB, AP, BP, AQ, BQ

    DZBW_CMP1 = 0_c_int
    IF (PA .EQ. PB) RETURN

    CALL C_F_POINTER(TRANSFER(PA, C_NULL_PTR), A)
    CALL C_F_POINTER(TRANSFER(PB, C_NULL_PTR), B)

    IF (A%W .LT. B%W) THEN
       DZBW_CMP1 = 1_c_int
    ELSE IF (A%W .GT. B%W) THEN
       DZBW_CMP1 = -1_c_int
    ELSE IF (A%W .EQ. B%W) THEN
       AB = ABS(A%B)
       BB = ABS(B%B)
       IF (AB .LT. BB) THEN
          DZBW_CMP1 = 2_c_int
       ELSE IF (AB .GT. BB) THEN
          DZBW_CMP1 = -2_c_int
       ELSE ! AB = BB
          AP = ABS(A%P)
          BP = ABS(B%P)
          IF (AP .LT. BP) THEN
             DZBW_CMP1 = 3_c_int
          ELSE IF (AP .GT. BP) THEN
             DZBW_CMP1 = -3_c_int
          ELSE ! AP = BP
             AQ = ABS(A%Q)
             BQ = ABS(B%Q)
             IF (AQ .LT. BQ) THEN
                DZBW_CMP1 = 4_c_int
             ELSE IF (AQ .GT. BQ) THEN
                DZBW_CMP1 = -4_c_int
             ELSE ! AQ = BQ
                DZBW_CMP1 = 0_c_int
             END IF
          END IF
       END IF
    ELSE ! NaN magnitude(s)
       IF (B%W .EQ. B%W) THEN
          DZBW_CMP1 = 5_c_int
       ELSE ! NaN(B%W)
          DZBW_CMP1 = -5_c_int
       END IF
    END IF
  END FUNCTION DZBW_CMP1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! sorting by type
  ! same type ==> sorting by magnitude
  !   same magnitude ==> sorting by subdiagonals (bands)
  !     outer bands first
  !       within a band, the lower elements first
  FUNCTION DZBW_CMP2(PA, PB)
    IMPLICIT NONE
    INTEGER(c_intptr_t), INTENT(IN), VALUE :: PA, PB
    INTEGER(c_int) :: DZBW_CMP2

    TYPE(DZBW), POINTER :: A, B
    INTEGER :: AB, BB, AP, BP, AQ, BQ

    DZBW_CMP2 = 0_c_int
    IF (PA .EQ. PB) RETURN

    CALL C_F_POINTER(TRANSFER(PA, C_NULL_PTR), A)
    CALL C_F_POINTER(TRANSFER(PB, C_NULL_PTR), B)

    IF (A%B .LT. 0) THEN
       IF (B%B .GE. 0) THEN
          DZBW_CMP2 = 1_c_int
          RETURN
       ELSE ! B%B < 0
          BB = -B%B
       END IF
       AB = -A%B
    ELSE ! A%B >= 0
       IF (B%B .LT. 0) THEN
          DZBW_CMP2 = -1_c_int
          RETURN
       ELSE ! B%B >= 0
          BB = B%B
       END IF
       AB = A%B
    END IF

    IF (A%W .LT. B%W) THEN
       DZBW_CMP2 = 2_c_int
    ELSE IF (A%W .GT. B%W) THEN
       DZBW_CMP2 = -2_c_int
    ELSE IF (A%W .EQ. B%W) THEN
       IF (AB .LT. BB) THEN
          DZBW_CMP2 = 3_c_int
       ELSE IF (AB .GT. BB) THEN
          DZBW_CMP2 = -3_c_int
       ELSE ! AB = BB
          AP = ABS(A%P)
          BP = ABS(B%P)
          IF (AP .LT. BP) THEN
             DZBW_CMP2 = 4_c_int
          ELSE IF (AP .GT. BP) THEN
             DZBW_CMP2 = -4_c_int
          ELSE ! AP = BP
             AQ = ABS(A%Q)
             BQ = ABS(B%Q)
             IF (AQ .LT. BQ) THEN
                DZBW_CMP2 = 5_c_int
             ELSE IF (AQ .GT. BQ) THEN
                DZBW_CMP2 = -5_c_int
             ELSE ! AQ = BQ
                DZBW_CMP2 = 0_c_int
             END IF
          END IF
       END IF
    ELSE ! NaN magnitude(s)
       IF (B%W .EQ. B%W) THEN
          DZBW_CMP2 = 6_c_int
       ELSE ! NaN(B%W)
          DZBW_CMP2 = -6_c_int
       END IF
    END IF
  END FUNCTION DZBW_CMP2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DZBW_SORT(NN, DZ, CMP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NN
    TYPE(DZBW), INTENT(INOUT), TARGET :: DZ(NN)
    PROCEDURE(VN_QSORT_CMP) :: CMP
    INTEGER, INTENT(OUT) :: INFO

    IF (NN .LT. 0) THEN
       INFO = -1
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN
    IF (NN .EQ. 0) RETURN

    INFO = GET_THREAD_NS()
    CALL VN_QSORT(C_LOC(DZ), INT(NN,c_size_t), C_SIZEOF(DZ(1)), C_FUNLOC(CMP))
    INFO = GET_THREAD_NS() - INFO
  END SUBROUTINE DZBW_SORT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DZBW_NCP(NN, DZ, N_2, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NN, N_2
    TYPE(DZBW), INTENT(IN) :: DZ(NN)
    INTEGER, INTENT(OUT) :: SL, STEP(N_2), INFO

    INTEGER :: I, J, K, AP, AQ, BP, BQ

    IF (NN .LT. 0) THEN
       INFO = -1
    ELSE IF (N_2 .LT. 0) THEN
       INFO = -3
    ELSE IF (N_2 .GT. NN) THEN
       INFO = -3
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    INFO = GET_THREAD_NS()
    SL = 0
    IF (NN .EQ. 0) GOTO 1
    IF (N_2 .EQ. 0) GOTO 1

    I = 1
    DO WHILE (SL .LT. N_2)
       SL = SL + 1
       STEP(SL) = I
       IF (SL .GE. N_2) EXIT

       AP = ABS(DZ(I)%P)
       AQ = ABS(DZ(I)%Q)
       K = NN + 1

       !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(J,BP,BQ) SHARED(I,NN,DZ,AP,AQ) REDUCTION(MIN:K)
       DO J = I+1, NN
          IF (DZ(J)%W .EQ. DZ(J)%W) THEN
             BP = ABS(DZ(J)%P)
             BQ = ABS(DZ(J)%Q)
             IF ((AP .EQ. BP) .OR. (AP .EQ. BQ) .OR. (AQ .EQ. BP) .OR. (AQ .EQ. BQ)) THEN
                DZ(J)%W = QUIET_NAN(J)
             ELSE ! not colliding
                K = MIN(K, J)
             END IF
          END IF
       END DO
       !$OMP END PARALLEL DO

       IF (K .GT. NN) EXIT
       I = K
    END DO

    !DIR$ VECTOR ALWAYS
    DO I = SL+1, N_2
       STEP(I) = 0
    END DO

1   INFO = GET_THREAD_NS() - INFO
  END SUBROUTINE DZBW_NCP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   SUBROUTINE DZBW_NCP(NN, DZ, N_2, SL, STEP, INFO)
!     IMPLICIT NONE
!     INTEGER, INTENT(IN) :: NN, N_2
!     TYPE(DZBW), INTENT(IN) :: DZ(NN)
!     INTEGER, INTENT(OUT) :: SL, STEP(N_2), INFO

!     INTEGER :: I, J, K, AP, AQ, BP, BQ
!     LOGICAL :: C

!     IF (NN .LT. 0) THEN
!        INFO = -1
!     ELSE IF (N_2 .LT. 0) THEN
!        INFO = -3
!     ELSE IF (N_2 .GT. NN) THEN
!        INFO = -3
!     ELSE ! all OK
!        INFO = 0
!     END IF
!     IF (INFO .NE. 0) RETURN

!     INFO = GET_THREAD_NS()
!     SL = 0
!     IF (NN .EQ. 0) GOTO 1
!     IF (N_2 .EQ. 0) GOTO 1

!     J = 1
!     DO I = 1, N_2
!        STEP(I) = J
!        J = J + 1
!        DO WHILE (J .LE. NN)
!           C = .NOT. (DZ(J)%W .EQ. DZ(J)%W)
!           IF (C) THEN
!              J = J + 1
!              CYCLE
!           END IF
!           AP = ABS(DZ(J)%P)
!           AQ = ABS(DZ(J)%Q)
!           DO K = I, 1, -1
!              BP = ABS(DZ(STEP(K))%P)
!              BQ = ABS(DZ(STEP(K))%Q)
!              IF ((AP .EQ. BP) .OR. (AP .EQ. BQ) .OR. (AQ .EQ. BP) .OR. (AQ .EQ. BQ)) THEN
!                 C = .TRUE.
!                 EXIT
!              END IF
!           END DO
!           IF (C) THEN
!              J = J + 1
!           ELSE
!              EXIT
!           END IF
!        END DO
!        SL = I
!        IF (J .GT. NN) EXIT
!     END DO

!     !DIR$ VECTOR ALWAYS
!     DO I = SL+1, N_2
!        STEP(I) = 0
!     END DO

! 1   INFO = GET_THREAD_NS() - INFO
!   END SUBROUTINE DZBW_NCP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE DTYPES
