MODULE ATYPES
  USE OMP_LIB
  USE TIMER
  USE VN_SORT_F
  IMPLICIT NONE

  ABSTRACT INTERFACE
     PURE INTEGER FUNCTION APQI(N, P, Q)
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: N, P, Q
     END FUNCTION APQI
  END INTERFACE

  ABSTRACT INTERFACE
     PURE SUBROUTINE ATRU(N, P, Q, NN, INFO)
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: N, NN
       INTEGER, INTENT(OUT) :: P(NN), Q(NN), INFO
     END SUBROUTINE ATRU
  END INTERFACE

  TYPE, BIND(C) :: AW
     REAL(KIND=DWP) :: W ! weight
     INTEGER :: P ! row
     INTEGER :: Q ! column
     INTEGER :: B ! band (Q - P) > 0
  END TYPE AW

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE AW_OUT(OU, HDR, NN, DZ, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: OU, NN, SL, STEP(SL)
    CHARACTER(LEN=*), INTENT(IN) :: HDR
    TYPE(AW), INTENT(IN) :: DZ(NN)
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

    INFO = GET_SYS_US()

    IF (LEN_TRIM(HDR) .GT. 0) THEN
       WRITE (OU,'(A)') TRIM(HDR)
    ELSE IF (SL .EQ. 0) THEN
       WRITE (OU,'(A)') '"J","W","P","Q","B"'
    ELSE ! SL > 0
       WRITE (OU,'(A)') '"I","J","W","P","Q","B"'
    END IF

    IF (SL .EQ. 0) THEN
       DO J = 1, NN
          WRITE (OU,'(I10,A,ES25.17E3,3(A,I10))') J, ',', DZ(J)%W, ',', DZ(J)%P, ',', DZ(J)%Q, ',', DZ(J)%B
       END DO
    ELSE ! SL > 0
       DO I = 1, SL
          J = STEP(I)
          WRITE (OU,'(2(I10,A),ES25.17E3,3(A,I10))') I, ',', J, ',', DZ(J)%W, ',', DZ(J)%P, ',', DZ(J)%Q, ',', DZ(J)%B
       END DO
    END IF

    INFO = GET_SYS_US() - INFO
  END SUBROUTINE AW_OUT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! sorting by magnitude
  !   same magnitude ==> sorting by subdiagonals (bands)
  !     outer bands first
  !       within a band, the lower elements first
#ifdef _GNU_SOURCE
  RECURSIVE INTEGER(KIND=c_int) FUNCTION AW_CMP1(PA, PB, CTX) BIND(C)
#else
  RECURSIVE INTEGER(KIND=c_int) FUNCTION AW_CMP1(CTX, PA, PB) BIND(C)
#endif
    IMPLICIT NONE
    TYPE(c_ptr), INTENT(IN), VALUE :: PA, PB, CTX

    TYPE(AW), POINTER :: A, B

    AW_CMP1 = 0_c_int
    CALL C_F_POINTER(PA, A)
    IF (.NOT. ASSOCIATED(A)) RETURN
    CALL C_F_POINTER(PB, B)
    IF (.NOT. ASSOCIATED(B)) RETURN
    IF (ASSOCIATED(A, B)) RETURN

    IF (A%W .LT. B%W) THEN
       AW_CMP1 = 1_c_int
    ELSE IF (A%W .GT. B%W) THEN
       AW_CMP1 = -1_c_int
    ELSE IF (A%W .EQ. B%W) THEN
       IF (A%B .LT. B%B) THEN
          AW_CMP1 = 2_c_int
       ELSE IF (A%B .GT. B%B) THEN
          AW_CMP1 = -2_c_int
       ELSE ! A%B = B%B
          IF (A%P .LT. B%P) THEN
             AW_CMP1 = 3_c_int
          ELSE IF (A%P .GT. B%P) THEN
             AW_CMP1 = -3_c_int
          ELSE ! A%P = B%P
             IF (A%Q .LT. B%Q) THEN
                AW_CMP1 = 4_c_int
             ELSE IF (A%Q .GT. B%Q) THEN
                AW_CMP1 = -4_c_int
             ELSE ! A%Q = B%Q
                AW_CMP1 = 0_c_int
             END IF
          END IF
       END IF
    ELSE ! NaN magnitude(s)
       IF (B%W .EQ. B%W) THEN
          AW_CMP1 = 5_c_int
       ELSE ! NaN(B%W)
          AW_CMP1 = -5_c_int
       END IF
    END IF
  END FUNCTION AW_CMP1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE AW_SRT1(NT, NN, NM, DZ, CMP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, NN, NM
    TYPE(AW), INTENT(INOUT), TARGET :: DZ(NM)
    PROCEDURE(VN_QSORT_CMP) :: CMP
    INTEGER, INTENT(OUT) :: INFO

    INTEGER :: T, TE, I, IE, J, JE, K, L, C, OE
    INTEGER :: EPA, EPT ! elements per array, thread
    TYPE(AW), POINTER, CONTIGUOUS :: A(:), B(:)

    IF (NT .LE. 1) THEN
       INFO = -1
    ELSE IF (NN .LT. 0) THEN
       INFO = -2
    ELSE IF (NM .LT. NN) THEN
       INFO = -3
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN
    IF (NN .EQ. 0) RETURN

    INFO = GET_SYS_US()

    EPA = NM / 2
    A => DZ(1:EPA)
    B => DZ(EPA+1:NM)
    EPT = EPA / NT

    ! virtual elements
    DO I = NN+1, EPA
       A(I)%W = QUIET_NAN(I)
    END DO

    ! Baudet-Stevenson odd-even sort with merge-splitting of the subarrays; see:
    ! Baudet and Stevenson, Optimal Sorting Algorithms for Parallel Computers,
    ! IEEE Transactions on Computers, C-27(1):84--87, Jan 1978.
    ! doi:10.1109/TC.1978.1674957

    !$OMP PARALLEL NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(I) SHARED(EPT,A)
    I = INT(OMP_GET_THREAD_NUM()) * EPT + 1
#ifdef _GNU_SOURCE
    CALL VN_QSORT(C_LOC(A(I)), INT(EPT,c_size_t), C_SIZEOF(A(I)), C_FUNLOC(CMP), C_NULL_PTR)
#else
    CALL VN_QSORT(C_LOC(A(I)), INT(EPT,c_size_t), C_SIZEOF(A(I)), C_NULL_PTR, C_FUNLOC(CMP))
#endif
    !$OMP END PARALLEL

    DO WHILE (.TRUE.)
       TE = 0
       DO OE = 0, 1
          L = 0
          !$OMP PARALLEL NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(T,I,IE,J,JE,K,C) SHARED(NT,EPT,OE,A,B) REDUCTION(+:L)
          T = INT(OMP_GET_THREAD_NUM())
          IF ((MOD(T, 2) .EQ. OE) .AND. ((T + 1) .LT. NT)) THEN
             ! merge with T + 1
             I = T * EPT + 1
             IE = I + (EPT - 1)
             J = IE + 1
             JE = J + (EPT - 1)
             K = I
             L = 0
             DO WHILE ((I .LE. IE) .AND. (J .LE. JE) .AND. (K .LE. JE))
#ifdef _GNU_SOURCE
                C = INT(CMP(C_LOC(A(I)), C_LOC(A(J)), C_NULL_PTR))
#else
                C = INT(CMP(C_NULL_PTR, C_LOC(A(I)), C_LOC(A(J))))
#endif
                IF (C .LE. 0) THEN
                   B(K) = A(I)
                   I = I + 1
                   K = K + 1
                ELSE ! C .GT. 0
                   B(K) = A(J)
                   J = J + 1
                   K = K + 1
                   L = L + 1
                END IF
             END DO
             DO WHILE ((I .LE. IE) .AND. (K .LE. JE))
                B(K) = A(I)
                I = I + 1
                K = K + 1
             END DO
             DO WHILE ((J .LE. JE) .AND. (K .LE. JE))
                B(K) = A(J)
                J = J + 1
                K = K + 1
             END DO
          END IF
          !$OMP END PARALLEL
          TE = TE + L
          !$OMP PARALLEL NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(I,IE,J) SHARED(EPT,A,B)
          I = INT(OMP_GET_THREAD_NUM()) * EPT + 1
          IE = I + (EPT - 1)
          DO J = I, IE
             A(J) = B(J)
          END DO
          !$OMP END PARALLEL
       END DO
       IF (TE .EQ. 0) EXIT
    END DO

    INFO = GET_SYS_US() - INFO
  END SUBROUTINE AW_SRT1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE AW_SRT2(NT, NN, NM, DZ, CMP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, NN, NM
    TYPE(AW), INTENT(INOUT), TARGET :: DZ(NM)
    PROCEDURE(VN_QSORT_CMP) :: CMP
    INTEGER, INTENT(OUT) :: INFO

    IF (NT .LE. 0) THEN
       INFO = -1
    ELSE IF (NN .LT. 0) THEN
       INFO = -2
    ELSE IF (NM .LT. NN) THEN
       INFO = -3
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN
    IF (NN .EQ. 0) RETURN

    INFO = GET_SYS_US()
#ifdef _GNU_SOURCE
    CALL VN_QSORT(C_LOC(DZ), INT(NN,c_size_t), C_SIZEOF(DZ(1)), C_FUNLOC(CMP), C_NULL_PTR)
#else
    CALL VN_QSORT(C_LOC(DZ), INT(NN,c_size_t), C_SIZEOF(DZ(1)), C_NULL_PTR, C_FUNLOC(CMP))
#endif
    INFO = GET_SYS_US() - INFO
  END SUBROUTINE AW_SRT2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RECURSIVE SUBROUTINE AW_NCP1(NT, N, NN, NM, DZ, N_2, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, N, NN, NM, N_2
    TYPE(AW), INTENT(INOUT) :: DZ(NM)
    INTEGER, INTENT(OUT) :: SL, STEP(N_2), INFO

    INTEGER :: I, J, K, AP, AQ, BP, BQ

    SL = 0
#ifndef NDEBUG
    DO I = 1, N_2
       STEP(I) = 0
    END DO
#endif

    IF (NT .LE. 0) THEN
       INFO = -1
    ELSE IF (N .LT. 0) THEN
       INFO = -2
    ELSE IF (NN .LT. 0) THEN
       INFO = -3
    ELSE IF (NM .LT. NN) THEN
       INFO = -4
    ELSE IF (N_2 .LT. 0) THEN
       INFO = -6
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    INFO = GET_SYS_US()
    IF (N .EQ. 0) GOTO 1
    IF (NN .EQ. 0) GOTO 1
    IF (N_2 .EQ. 0) GOTO 1

    I = 1
    DO WHILE (SL .LT. N_2)
       SL = SL + 1
       STEP(SL) = I
       IF (SL .GE. N_2) EXIT

       AP = DZ(I)%P
       AQ = DZ(I)%Q
       K = NN + 1

       !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(J,BP,BQ) SHARED(I,NN,DZ,AP,AQ) REDUCTION(MIN:K)
       DO J = I+1, NN
          IF (DZ(J)%W .EQ. DZ(J)%W) THEN
             BP = DZ(J)%P
             BQ = DZ(J)%Q
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

1   INFO = GET_SYS_US() - INFO
  END SUBROUTINE AW_NCP1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RECURSIVE SUBROUTINE AW_NCP2(NT, N, NN, NM, DZ, N_2, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, N, NN, NM, N_2
    TYPE(AW), INTENT(INOUT) :: DZ(NM)
    INTEGER, INTENT(OUT) :: SL, STEP(N_2), INFO

    INTEGER :: I, J, K, AP, AQ, BP, BQ
    LOGICAL :: C

    SL = 0
#ifndef NDEBUG
    DO I = 1, N_2
       STEP(I) = 0
    END DO
#endif

    IF (NT .LE. 0) THEN
       INFO = -1
    ELSE IF (N .LT. 0) THEN
       INFO = -2
    ELSE IF (NN .LT. 0) THEN
       INFO = -3
    ELSE IF (NM .LT. NN) THEN
       INFO = -4
    ELSE IF (N_2 .LT. 0) THEN
       INFO = -6
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    INFO = GET_SYS_US()
    IF (N .EQ. 0) GOTO 2
    IF (NN .EQ. 0) GOTO 2
    IF (N_2 .EQ. 0) GOTO 2

    J = 1
    DO I = 1, N_2
       STEP(I) = J
       J = J + 1
       DO WHILE (J .LE. NN)
          C = .NOT. (DZ(J)%W .EQ. DZ(J)%W)
          IF (C) THEN
             J = J + 1
             CYCLE
          END IF
          AP = DZ(J)%P
          AQ = DZ(J)%Q
          DO K = I, 1, -1
             BP = DZ(STEP(K))%P
             BQ = DZ(STEP(K))%Q
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
       SL = I
       IF (J .GT. NN) EXIT
    END DO

2   INFO = GET_SYS_US() - INFO
  END SUBROUTINE AW_NCP2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  RECURSIVE SUBROUTINE AW_NCP3(NT, NN, N, NM, DZ, N_2, SL, STEP, INFO)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, N, NN, NM, N_2
    TYPE(AW), INTENT(INOUT) :: DZ(NM)
    INTEGER, INTENT(OUT) :: SL, STEP(N_2), INFO

    INTEGER, TARGET :: STP(N_2)
    INTEGER :: I, J, K, L, M, AP, AQ, BP, BQ
    LOGICAL :: C
    REAL(KIND=DWP) :: W1, WL

    SL = 0
#ifndef NDEBUG
    DO I = 1, N_2
       STEP(I) = 0
    END DO
#endif

    IF (NT .LE. 0) THEN
       INFO = -1
    ELSE IF (N .LT. 0) THEN
       INFO = -2
    ELSE IF (NN .LT. 0) THEN
       INFO = -3
    ELSE IF (NM .LT. NN) THEN
       INFO = -4
    ELSE IF (N_2 .LT. 0) THEN
       INFO = -6
    ELSE ! all OK
       INFO = 0
    END IF
    IF (INFO .NE. 0) RETURN

    INFO = GET_SYS_US()
    IF (N .EQ. 0) GOTO 3
    IF (NN .EQ. 0) GOTO 3
    IF (N_2 .EQ. 0) GOTO 3

    J = 1
    DO I = 1, N_2
       STEP(I) = J
       J = J + 1
       DO WHILE (J .LE. NN)
          C = .NOT. (DZ(J)%W .EQ. DZ(J)%W)
          IF (C) THEN
             J = J + 1
             CYCLE
          END IF
          AP = DZ(J)%P
          AQ = DZ(J)%Q
          DO K = I, 1, -1
             BP = DZ(STEP(K))%P
             BQ = DZ(STEP(K))%Q
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
       SL = I
       IF (J .GT. NN) EXIT
    END DO
    IF (SL .EQ. 0) RETURN
    W1 = D_ZERO
    DO I = SL, 1, -1
       W1 = W1 + DZ(STEP(I))%W
    END DO

#ifdef OLD_OMP
    !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(I,J,K,L,M,AP,AQ,BP,BQ,C,WL,STP) SHARED(NN,N_2,SL,DZ,STEP,W1) &
    !$OMP& SCHEDULE(DYNAMIC,1)
#else
    !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) PRIVATE(I,J,K,L,M,AP,AQ,BP,BQ,C,WL,STP) SHARED(NN,N_2,SL,DZ,STEP,W1) &
    !$OMP& SCHEDULE(NONMONOTONIC:DYNAMIC,1)
#endif
    DO L = 2, NN
       M = 0
       J = L
       DO I = 1, N_2
          STP(I) = J
          J = J + 1
          DO WHILE (J .LE. NN)
             C = .NOT. (DZ(J)%W .EQ. DZ(J)%W)
             IF (C) THEN
                J = J + 1
                CYCLE
             END IF
             AP = DZ(J)%P
             AQ = DZ(J)%Q
             DO K = I, 1, -1
                BP = DZ(STP(K))%P
                BQ = DZ(STP(K))%Q
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
          M = I
          IF (J .GT. NN) EXIT
       END DO
       WL = D_ZERO
       DO I = M, 1, -1
          WL = WL + DZ(STP(I))%W
       END DO
       !$OMP CRITICAL
       IF (.NOT. (WL .LE. W1)) THEN
          SL = M
          DO I = 1, SL
             STEP(I) = STP(I)
          END DO
          W1 = WL
       END IF
       !$OMP END CRITICAL
    END DO
    !$OMP END PARALLEL DO

3   INFO = GET_SYS_US() - INFO
  END SUBROUTINE AW_NCP3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE ATYPES
