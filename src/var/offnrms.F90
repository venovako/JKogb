!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  REAL(KIND=DWP) FUNCTION DSTEP_OFFNORM1(NT, M, N, A, LDA)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, M, N, LDA
    REAL(KIND=DWP), INTENT(IN) :: A(LDA,N)

    REAL(KIND=DWP) :: S
    INTEGER :: I, J

    IF (NT .LE. 0) THEN
       DSTEP_OFFNORM1 = -1
    ELSE IF (M .LT. 0) THEN
       DSTEP_OFFNORM1 = -2
    ELSE IF (M .GT. LDA) THEN
       DSTEP_OFFNORM1 = -2
    ELSE IF (N .LT. 0) THEN
       DSTEP_OFFNORM1 = -3
    ELSE IF (N .GT. M) THEN
       DSTEP_OFFNORM1 = -3
    ELSE IF (LDA .LT. 0) THEN
       DSTEP_OFFNORM1 = -5
    ELSE ! all OK
       DSTEP_OFFNORM1 = 0
    END IF
    IF (DSTEP_OFFNORM1 .NE. D_ZERO) RETURN

    S = D_ZERO
    !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) SHARED(M,N,A) PRIVATE(I,J) REDUCTION(+:S)
    DO J = 1, N
       !DIR$ VECTOR ALWAYS
       DO I = 1, J-1
          S = S + ABS(A(I,J))
       END DO
       !DIR$ VECTOR ALWAYS
       DO I = J+1, M
          S = S + ABS(A(I,J))
       END DO
    END DO
    !$OMP END PARALLEL DO
    DSTEP_OFFNORM1 = S
  END FUNCTION DSTEP_OFFNORM1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  REAL(KIND=DWP) FUNCTION DSTEP_OFFNORMI(NT, M, N, A, LDA)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, M, N, LDA
    REAL(KIND=DWP), INTENT(IN) :: A(LDA,N)

    REAL(KIND=DWP) :: S
    INTEGER :: I, J

    IF (NT .LE. 0) THEN
       DSTEP_OFFNORMI = -1
    ELSE IF (M .LT. 0) THEN
       DSTEP_OFFNORMI = -2
    ELSE IF (M .GT. LDA) THEN
       DSTEP_OFFNORMI = -2
    ELSE IF (N .LT. 0) THEN
       DSTEP_OFFNORMI = -3
    ELSE IF (N .GT. M) THEN
       DSTEP_OFFNORMI = -3
    ELSE IF (LDA .LT. 0) THEN
       DSTEP_OFFNORMI = -5
    ELSE ! all OK
       DSTEP_OFFNORMI = 0
    END IF
    IF (DSTEP_OFFNORMI .NE. D_ZERO) RETURN

    S = D_ZERO
    !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) SHARED(M,N,A) PRIVATE(I,J) REDUCTION(MAX:S)
    DO J = 1, N
       !DIR$ VECTOR ALWAYS
       DO I = 1, J-1
          S = MAX(ABS(A(I,J)) + ABS(A(J,I)), S)
       END DO
    END DO
    !$OMP END PARALLEL DO
    DSTEP_OFFNORMI = S
  END FUNCTION DSTEP_OFFNORMI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  REAL(KIND=DWP) FUNCTION DSTEP_RELNORM1(NT, M, N, A, LDA)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, M, N, LDA
    REAL(KIND=DWP), INTENT(IN) :: A(LDA,N)

    REAL(KIND=DWP) :: OFFNORM, MINDIAG
    INTEGER :: I

    OFFNORM = DSTEP_OFFNORM1(NT, M, N, A, LDA)
    IF ((OFFNORM .LE. D_ZERO) .OR. .NOT. (OFFNORM .LE. HUGE(OFFNORM))) THEN
       DSTEP_RELNORM1 = OFFNORM
       RETURN
    END IF

    MINDIAG = HUGE(MINDIAG)
    !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) SHARED(N,A) PRIVATE(I) REDUCTION(MIN:MINDIAG)
    DO I = 1, N
       MINDIAG = MIN(ABS(A(I,I)), MINDIAG)
    END DO
    !$OMP END PARALLEL DO
    DSTEP_RELNORM1 = OFFNORM / MINDIAG
  END FUNCTION DSTEP_RELNORM1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  REAL(KIND=DWP) FUNCTION DSTEP_RELNORMI(NT, M, N, A, LDA)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, M, N, LDA
    REAL(KIND=DWP), INTENT(IN) :: A(LDA,N)

    REAL(KIND=DWP) :: OFFNORM, MINDIAG
    INTEGER :: I

    OFFNORM = DSTEP_OFFNORMI(NT, M, N, A, LDA)
    IF ((OFFNORM .LE. D_ZERO) .OR. .NOT. (OFFNORM .LE. HUGE(OFFNORM))) THEN
       DSTEP_RELNORMI = OFFNORM
       RETURN
    END IF

    MINDIAG = HUGE(MINDIAG)
    !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) SHARED(N,A) PRIVATE(I) REDUCTION(MIN:MINDIAG)
    DO I = 1, N
       MINDIAG = MIN(ABS(A(I,I)), MINDIAG)
    END DO
    !$OMP END PARALLEL DO
    DSTEP_RELNORMI = OFFNORM / MINDIAG
  END FUNCTION DSTEP_RELNORMI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! IF (Q .EQ. 0) THEN
       !    TW = DSTEP_RELNORMI(NT, N, N, A, LDA)
       !    IF (TW .LT. EPSILON(TW)) THEN
       !       I = -1
       !    ELSE ! not converged in this sense
       !       I = 1
       !    END IF
       ! ELSE ! diag(A) has changed
       !    I = 0
       ! END IF

       ! IF (I .EQ. -1) INFO = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  REAL(KIND=DWP) FUNCTION ZSTEP_OFFNORM1(NT, M, N, A, LDA)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, M, N, LDA
    COMPLEX(KIND=DWP), INTENT(IN) :: A(LDA,N)

    REAL(KIND=DWP) :: S
    INTEGER :: I, J

    IF (NT .LE. 0) THEN
       ZSTEP_OFFNORM1 = -1
    ELSE IF (M .LT. 0) THEN
       ZSTEP_OFFNORM1 = -2
    ELSE IF (M .GT. LDA) THEN
       ZSTEP_OFFNORM1 = -2
    ELSE IF (N .LT. 0) THEN
       ZSTEP_OFFNORM1 = -3
    ELSE IF (N .GT. M) THEN
       ZSTEP_OFFNORM1 = -3
    ELSE IF (LDA .LT. 0) THEN
       ZSTEP_OFFNORM1 = -5
    ELSE ! all OK
       ZSTEP_OFFNORM1 = 0
    END IF
    IF (ZSTEP_OFFNORM1 .NE. D_ZERO) RETURN

    S = D_ZERO
    !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) SHARED(M,N,A) PRIVATE(I,J) REDUCTION(+:S)
    DO J = 1, N
       !DIR$ VECTOR ALWAYS
       DO I = 1, J-1
          S = S + ABS(A(I,J))
       END DO
       !DIR$ VECTOR ALWAYS
       DO I = J+1, M
          S = S + ABS(A(I,J))
       END DO
    END DO
    !$OMP END PARALLEL DO
    ZSTEP_OFFNORM1 = S
  END FUNCTION ZSTEP_OFFNORM1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  REAL(KIND=DWP) FUNCTION ZSTEP_RELNORM1(NT, M, N, A, LDA)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT, M, N, LDA
    COMPLEX(KIND=DWP), INTENT(IN) :: A(LDA,N)

    REAL(KIND=DWP) :: OFFNORM, MINDIAG
    INTEGER :: I

    OFFNORM = ZSTEP_OFFNORM1(NT, M, N, A, LDA)
    IF ((OFFNORM .LE. D_ZERO) .OR. .NOT. (OFFNORM .LE. HUGE(OFFNORM))) THEN
       ZSTEP_RELNORM1 = OFFNORM
       RETURN
    END IF

    MINDIAG = HUGE(MINDIAG)
    !$OMP PARALLEL DO NUM_THREADS(NT) DEFAULT(NONE) SHARED(N,A) PRIVATE(I) REDUCTION(MIN:MINDIAG)
    DO I = 1, N
       MINDIAG = MIN(ABS(A(I,I)), MINDIAG)
    END DO
    !$OMP END PARALLEL DO
    ZSTEP_RELNORM1 = OFFNORM / MINDIAG
  END FUNCTION ZSTEP_RELNORM1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! IF (Q .EQ. 0) THEN
       !    TW = ZSTEP_RELNORM1(NT, N, N, A, LDA)
       !    IF (TW .LT. EPSILON(TW)) THEN
       !       I = -1
       !    ELSE ! not converged in this sense
       !       I = 1
       !    END IF
       ! ELSE ! diag(A) has changed
       !    I = 0
       ! END IF

       ! IF (I .EQ. -1) INFO = 0
