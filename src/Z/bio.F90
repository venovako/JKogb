!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ZOPEN_YJ_RO(FN, M, N, SZ, FD, INFO)
  IMPLICIT NONE
  CHARACTER(LEN=*,KIND=c_char), INTENT(IN) :: FN
  INTEGER, INTENT(IN) :: M, N
  INTEGER, INTENT(OUT) :: SZ(2), FD(2), INFO

  INTEGER :: EXPTSZ(2), DIFFSZ(2)

  SZ = -1
  FD = -1

  IF (M .LT. 0) THEN
     INFO = -2
  ELSE IF (N .LT. 0) THEN
     INFO = -3
  ELSE
     INFO = 0
  END IF
  IF (INFO .NE. 0) RETURN

  EXPTSZ(1) = M * N * C_SIZEOF(Z_ZERO)
  EXPTSZ(2) = M     * C_SIZEOF(0)

  DIFFSZ = 0

  CALL BOPEN_RO((TRIM(FN)//c_char_'.Y'), SZ(1), FD(1))
  IF (FD(1) .LT. 0) THEN
     INFO = 1
     GOTO 1
  END IF
  DIFFSZ(1) = SZ(1) - EXPTSZ(1)
  IF (DIFFSZ(1) .NE. 0) THEN
     INFO = 1
     GOTO 1
  END IF
 
  CALL BOPEN_RO((TRIM(FN)//c_char_'.J'), SZ(2), FD(2))
  IF (FD(2) .LT. 0) THEN
     INFO = 2
     GOTO 1
  END IF
  DIFFSZ(2) = SZ(2) - EXPTSZ(2)
  IF (DIFFSZ(2) .NE. 0) THEN
     INFO = 2
     GOTO 1
  END IF

  RETURN

1 SZ = DIFFSZ
END SUBROUTINE ZOPEN_YJ_RO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ZREAD_YJ(FD, Y, J, M, N, SZ, INFO)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: FD(2), M, N
  COMPLEX(KIND=DWP), INTENT(OUT), TARGET :: Y(M,N)
  INTEGER, INTENT(OUT), TARGET :: J(M)
  INTEGER, INTENT(OUT) :: SZ(2), INFO

  INFO = 0

  SZ(1) = BREAD(FD(1), C_LOC(Y), C_SIZEOF(Y), 0)
  IF (SZ(1) .NE. C_SIZEOF(Y)) INFO = 1

  SZ(2) = BREAD(FD(2), C_LOC(J), C_SIZEOF(J), 0)
  IF (SZ(2) .NE. C_SIZEOF(J)) INFO = 2
END SUBROUTINE ZREAD_YJ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ZOPEN_UZS_RW(FN, M, N, SZ, FD, INFO)
  IMPLICIT NONE
  CHARACTER(LEN=*,KIND=c_char), INTENT(IN) :: FN
  INTEGER, INTENT(IN) :: M, N
  INTEGER, INTENT(OUT) :: SZ(3), FD(3), INFO

  SZ = -1
  FD = -1

  IF (M .LT. 0) THEN
     INFO = -2
  ELSE IF (N .LT. 0) THEN
     INFO = -3
  ELSE
     INFO = 0
  END IF
  IF (INFO .NE. 0) RETURN

  SZ(1) = M * N * C_SIZEOF(Z_ZERO) ! YU
  SZ(2) = N * N * C_SIZEOF(Z_ZERO) ! ZZ
  SZ(3) =     N * C_SIZEOF(D_ZERO) ! SY

  CALL BOPEN_RW((TRIM(FN)//c_char_'.YU'), SZ(1), FD(1))
  IF (FD(1) .LT. 0) THEN
     INFO = 1
     RETURN
  END IF

  CALL BOPEN_RW((TRIM(FN)//c_char_'.ZZ'), SZ(2), FD(2))
  IF (FD(2) .LT. 0) THEN
     INFO = 2
     RETURN
  END IF

  CALL BOPEN_RW((TRIM(FN)//c_char_'.SY'), SZ(3), FD(3))
  IF (FD(3) .LT. 0) THEN
     INFO = 3
     RETURN
  END IF
END SUBROUTINE ZOPEN_UZS_RW

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ZWRITE_UZS(FD, YU, Z, SY, M, N, INFO)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: FD(3), M, N
  COMPLEX(KIND=DWP), INTENT(IN), TARGET :: YU(M,N), Z(N,N)
  REAL(KIND=DWP), INTENT(IN), TARGET :: SY(N)
  INTEGER, INTENT(OUT) :: INFO

  INTEGER(KIND=c_size_t) :: SZ(3), S

  SZ(1) = M * N * C_SIZEOF(Z_ZERO) ! YU
  SZ(2) = N * N * C_SIZEOF(Z_ZERO) !  Z
  SZ(3) =     N * C_SIZEOF(D_ZERO) ! SY

  INFO = 0

  S = BWRITE(FD(1), C_LOC(YU), C_SIZEOF(YU), 0)
  IF (S .NE. SZ(1)) INFO = 1
  S = BWRITE(FD(2), C_LOC(Z), C_SIZEOF(Z), 0)
  IF (S .NE. SZ(2)) INFO = 2
  S = BWRITE(FD(3), C_LOC(SY), C_SIZEOF(SY), 0)
  IF (S .NE. SZ(3)) INFO = 3
END SUBROUTINE ZWRITE_UZS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
