PROGRAM DTEST
  USE DTRANSF
  IMPLICIT NONE

  REAL(KIND=DWP) :: A11, A21, A12, A22, CU, SU, CZ, SZ, A(2,2), U(2,2), Z(2,2)
  INTEGER :: J1, J2
  LOGICAL :: LEFT, LI, RI, NC

  WRITE (*,'(A)',ADVANCE='NO') 'A11='
  READ (*,*) A11
  WRITE (*,'(A)',ADVANCE='NO') 'A21='
  READ (*,*) A21
  WRITE (*,'(A)',ADVANCE='NO') 'A12='
  READ (*,*) A12
  WRITE (*,'(A)',ADVANCE='NO') 'A22='
  READ (*,*) A22
  WRITE (*,'(A)',ADVANCE='NO') 'J1='
  READ (*,*) J1
  WRITE (*,'(A)',ADVANCE='NO') 'J2='
  READ (*,*) J2
  WRITE (*,'(A)',ADVANCE='NO') 'LEFT='
  READ (*,*) LEFT

  CALL DTRANSFA(A11, A21, A12, A22, J1, J2, CU, SU, CZ, SZ, LEFT, A, U, Z, LI, RI, NC)

  WRITE (*,*) 'LI=', LI
  WRITE (*,*) 'RI=', RI
  WRITE (*,*) 'NC=', NC

  WRITE (*,*) 'CU=', CU
  WRITE (*,*) 'SU=', SU
  WRITE (*,*) 'CZ=', CZ
  WRITE (*,*) 'SZ=', SZ
  
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
