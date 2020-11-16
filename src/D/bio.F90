!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE DREAD_YJ(FD, Y, J, M, N, INFO)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: FD(2), M, N
  REAL(KIND=DWP), INTENT(OUT) :: Y(M,N)
  INTEGER, INTENT(OUT) :: J(M), INFO

  CALL BIO_READ_D2(FD(1), M, N, Y, M, INFO)
  IF (INFO .NE. 0) THEN
     INFO = 1
     RETURN
  END IF

  CALL BIO_READ_I1(FD(2), M, J, INFO)
  IF (INFO .NE. 0) THEN
     INFO = 2
     RETURN
  END IF
END SUBROUTINE DREAD_YJ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE DWRITE_UZS(FD, YU, ZZ, SY, M, N, INFO)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: FD(3), M, N
  REAL(KIND=DWP), INTENT(IN) :: YU(M,N), ZZ(N,N), SY(N)
  INTEGER, INTENT(OUT) :: INFO

  CALL BIO_WRITE_D2(FD(1), M, N, YU, M, INFO)
  IF (INFO .NE. 0) THEN
     INFO = 1
     RETURN
  END IF

  CALL BIO_WRITE_D2(FD(2), N, N, ZZ, N, INFO)
  IF (INFO .NE. 0) THEN
     INFO = 2
     RETURN
  END IF

  CALL BIO_WRITE_D1(FD(3), N, SY, INFO)
  IF (INFO .NE. 0) THEN
     INFO = 3
     RETURN
  END IF
END SUBROUTINE DWRITE_UZS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
