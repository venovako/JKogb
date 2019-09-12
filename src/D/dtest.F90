PROGRAM DTEST
  USE DTRANSF
  IMPLICIT NONE

  REAL(KIND=DWP) :: A(2,2), U(2,2), Z(2,2)
  INTEGER :: J(2)

  WRITE (*,'(A)',ADVANCE='NO') 'A11='
  READ (*,*) A(1,1)
  WRITE (*,'(A)',ADVANCE='NO') 'A21='
  READ (*,*) A(2,1)
  WRITE (*,'(A)',ADVANCE='NO') 'A12='
  READ (*,*) A(1,2)
  WRITE (*,'(A)',ADVANCE='NO') 'A22='
  READ (*,*) A(2,2)
  WRITE (*,'(A)',ADVANCE='NO') 'J1='
  READ (*,*) J(1)
  WRITE (*,'(A)',ADVANCE='NO') 'J2='
  READ (*,*) J(2)

  CALL DTRANSFA(A, J, U, Z)

  WRITE (*,*) 'A11=', A(1,1)
  WRITE (*,*) 'A21=', A(2,1)
  WRITE (*,*) 'A12=', A(1,2)
  WRITE (*,*) 'A22=', A(2,2)

  WRITE (*,*) 'U11=', U(1,1)
  WRITE (*,*) 'U21=', U(2,1)
  WRITE (*,*) 'U12=', U(1,2)
  WRITE (*,*) 'U22=', U(2,2)

  WRITE (*,*) 'Z11=', Z(1,1)
  WRITE (*,*) 'Z21=', Z(2,1)
  WRITE (*,*) 'Z12=', Z(1,2)
  WRITE (*,*) 'Z22=', Z(2,2)
END PROGRAM DTEST
