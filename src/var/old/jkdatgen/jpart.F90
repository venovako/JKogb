SUBROUTINE JPART(NROW, NCOLR, G, LDG, JVEC, NPLUS, IPL, INVP)

  !     Purpose
  !     =======
  !
  !     Transforms the G*J*G^T factorization into G_1*J_part*G_1^T
  !     factorization, with J_part partitioned as  J_part = ( I, -I ).
  !     Reorders the columns of G and the elements of JVEC.
  !     NPLUS is the number of elements in JVEC equal to 1.

  IMPLICIT NONE
#include "wp.F"

  INTEGER, INTENT(IN) :: NROW, NCOLR, LDG
  REAL(KIND=WP), DIMENSION(LDG, NROW), INTENT(INOUT) :: G
  INTEGER, DIMENSION(NROW), INTENT(INOUT) :: JVEC
  INTEGER, INTENT(OUT) :: NPLUS, IPL(NROW), INVP(NROW)

  INTEGER :: I, IPLUS, IMINUS, IP, JTEMP

  EXTERNAL :: XSWAP

  !     Count columns with JVEC( I ) = 1.
  NPLUS = 0
  DO I = 1, NCOLR
     IF ( JVEC( I ) .EQ. 1 ) NPLUS = NPLUS + 1
  END DO

  !     Set permutation IPL, where IPL( I ) holds the current place
  !     of the final I-th column.
  !     The following algorithm preserves the relative order of columns
  !     with the same sign in JVEC.
  IPLUS = 0
  IMINUS = NPLUS
  DO I = 1, NCOLR
     IF ( JVEC( I ) .EQ. 1 ) THEN
        IPLUS = IPLUS + 1
        IPL( IPLUS ) = I
     ELSE
        IMINUS = IMINUS + 1
        IPL( IMINUS ) = I
     END IF
  END DO
  DO I = NCOLR + 1, NROW
     IPL( I ) = I
  END DO

  !     Invert the permutation IPL and store it in INVP.
  DO I = 1, NROW
     INVP( IPL( I ) ) = I
  END DO

  !     Early return - all JVEC( I ) have the same sign.
  IF ( ( NPLUS .EQ. 0 ) .OR. ( NPLUS .EQ. NCOLR ) ) GOTO 1

  DO I = 1, NCOLR
     !     Swap columns G( I ) and G( IPL( I ) ).
     !     Also swap the corresponding elements in JVEC.
     IF ( IPL( I ) .NE. I ) THEN
        IP = IPL( I )

        CALL XSWAP( NROW, G( 1, I ), 1, G( 1, IP ), 1 )
        JTEMP = JVEC( I )
        JVEC( I ) = JVEC( IP )
        JVEC( IP ) = JTEMP

        INVP( IP ) = INVP( I )
        IPL( INVP( I ) ) = IP
     END IF
     !     Not necessary to set:
     !     IPL( I ) = I
     !     INVP( I ) = I
  END DO

1 DO I = NCOLR + 1, NROW
     JVEC( I ) = 0
  END DO

END SUBROUTINE JPART
