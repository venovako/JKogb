!  Pointwise J-Kogbetliantz
MODULE JK1
  IMPLICIT NONE

  INCLUDE 'params.f90'
  INCLUDE 'c_ifcs.f90'

CONTAINS

  INCLUDE 'dkogul.f90'
  INCLUDE 'dkogt2.f90'
  INCLUDE 'dkogh2.f90'
  INCLUDE 'jkstep.f90'

  SUBROUTINE DJKSVD(UPLO, M, N, G, LDG, U, LDU, V, LDV, JVEC, NPLUS, TOL, NSWEEP, IFLAGS, INFO)
    IMPLICIT NONE

    !   nsweep .eq. 0: iterate until convergence

    CHARACTER(LEN=1), INTENT(IN) :: UPLO
    INTEGER, INTENT(IN) :: M, N, LDG, LDU, LDV, JVEC(N), NPLUS, NSWEEP, IFLAGS
    DOUBLE PRECISION, INTENT(IN) :: TOL
    DOUBLE PRECISION, INTENT(INOUT) :: G(LDG,N), U(LDU,M), V(LDV,N)
    INTEGER, INTENT(OUT) :: INFO(2)

    INTEGER :: MAXCYC, SWEEP, SWHALF, STEP
    INTEGER :: I, J, STATUS
    LOGICAL :: UPPER, ROWCYC, HYP

    EXTERNAL :: DSCAL, DSWAP
    LOGICAL, EXTERNAL :: LSAME

    IF (LSAME(UPLO, 'U')) THEN
       UPPER = .TRUE.
    ELSE IF (LSAME(UPLO, 'L')) THEN
       UPPER = .FALSE.
    ELSE
       INFO(1) = -1
       INFO(2) = ICHAR(UPLO)
       RETURN
    END IF

    IF (M .LT. 0) THEN
       INFO(1) = -2
       INFO(2) = M
       RETURN
    END IF

    IF ((N .LT. 0) .OR. (N .GT. M)) THEN
       INFO(1) = -3
       INFO(2) = N
       RETURN
    END IF
    
    IF (LDG .LT. M) THEN
       INFO(1) = -5
       INFO(2) = LDG
       RETURN
    END IF
    
    IF (LDU .LT. M) THEN
       INFO(1) = -7
       INFO(2) = LDU
       RETURN
    END IF

    IF (LDV .LT. N) THEN
       INFO(1) = -9
       INFO(2) = LDV
       RETURN
    END IF

    IF ((NPLUS .LT. 0) .OR. (NPLUS .GT. N)) THEN
       INFO(1) = -11
       INFO(2) = NPLUS
       RETURN
    END IF

    IF (TOL .LT. ZERO) THEN
       INFO(1) = -12
       INFO(2) = -1
       RETURN
    END IF

    IF (NSWEEP .LT. 0) THEN
       INFO(1) = -13
       INFO(2) = NSWEEP
       RETURN
    END IF

    IF ((IFLAGS .LT. 0) .OR. (IFLAGS .GE. IFLAG_MAXFLG)) THEN
       INFO(1) = -14
       INFO(2) = IFLAGS
       RETURN
    END IF

    INFO = 0

    IF (M .EQ. 0) RETURN

    ROWCYC = (IAND(IFLAGS, IFLAG_ROWCYC) .NE. 0)
    IF (NSWEEP .EQ. 0) THEN
       MAXCYC = HUGE(MAXCYC)
    ELSE
       MAXCYC = NSWEEP
    END IF

    IF (ROWCYC) THEN
       DO SWEEP = 1, MAXCYC
          DO SWHALF = 0, 1
             STEP = 0
             STATUS = 0

             WRITE (*,9) 'Sweep: ', SWEEP, ' swhalf: ', SWHALF

             DO I = 1, N-1
                DO J = I+1, N
                   CALL JKSTEP(UPPER, I, J, M, N, G, LDG, U, LDU, V, LDV, JVEC, NPLUS, TOL, SWEEP, SWHALF, STEP, HYP, INFO)
                   IF (INFO(1) .NE. 0) RETURN
                   STATUS = MAX(STATUS, INFO(2))
                END DO
             END DO
             IF (STATUS .EQ. 0) GOTO 1
             UPPER = .NOT. UPPER
          END DO
       END DO
    ELSE
       DO SWEEP = 1, MAXCYC
          DO SWHALF = 0, 1
             STEP = 0
             STATUS = 0
             
             WRITE (*,9) 'Sweep: ', SWEEP, ' swhalf: ', SWHALF

             DO J = 2, N
                DO I = 1, J-1
                   CALL JKSTEP(UPPER, I, J, M, N, G, LDG, U, LDU, V, LDV, JVEC, NPLUS, TOL, SWEEP, SWHALF, STEP, HYP, INFO)
                   IF (INFO(1) .NE. 0) RETURN
                   STATUS = MAX(STATUS, INFO(2))
                END DO
             END DO
             IF (STATUS .EQ. 0) GOTO 1
             UPPER = .NOT. UPPER
          END DO
       END DO
    END IF

1   IF (SWEEP .LE. MAXCYC) THEN
       INFO(2) = 2 * SWEEP + SWHALF - 1
    ELSE
       INFO(2) = -MAXCYC
    END IF

    IF (IAND(IFLAGS, IFLAG_PPROCU) .NE. 0) THEN
       DO J = 1, M-1
          CALL DSWAP(M-J, U(J+1,J), 1, U(J, J+1), LDU)
       END DO
    END IF

    IF (IAND(IFLAGS, IFLAG_PPROCV) .NE. 0) THEN
       DO J = 1, N
          IF (JVEC(J) .NE. 1) CALL DSCAL(N, DBLE(JVEC(J)), V(1,J), 1)
       END DO
       DO I = 1, N
          IF (JVEC(I) .NE. 1) CALL DSCAL(N, DBLE(JVEC(I)), V(I,1), LDV)
       END DO
       DO J = 1, N-1
          CALL DSWAP(N-J, V(J+1,J), 1, V(J,J+1), LDV)
       END DO
    END IF

9   FORMAT (A,I2.2,A,I1.1)
  END SUBROUTINE DJKSVD
END MODULE JK1
