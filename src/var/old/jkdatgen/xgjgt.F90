SUBROUTINE XGJGT( INFO, N, NRANK, G, LDA, IPIV, JVEC, INFOG )

  !     Purpose
  !     =======
  !
  !     GJGT computes the indefinite factorization from XSYTRF.
  !     Returns G in the LOWER triangle of G.

  IMPLICIT NONE
#include "wp.F"

  INTEGER, INTENT(IN) :: INFO
  INTEGER, INTENT(IN) :: N
  INTEGER, INTENT(INOUT) :: NRANK
  REAL(KIND=WP), DIMENSION(LDA, N), INTENT(INOUT) :: G
  INTEGER, INTENT(IN) :: LDA
  INTEGER, DIMENSION(N), INTENT(INOUT) :: IPIV
  INTEGER, DIMENSION(N), INTENT(OUT) :: JVEC
  INTEGER, INTENT(OUT) :: INFOG

  REAL(KIND=WP), PARAMETER :: ZERO = +0.0E+0_WP

  INTEGER :: NR2, K, KP
  REAL(KIND=WP) :: AA, BB, CC, CS1, SN1, D, E

  EXTERNAL :: XLAEV2, XLASET, XROT, XSCAL, XSWAP

  INFOG = 0
  IF ( INFO .LT. 0 ) RETURN

  IF ( NRANK .LT. 0 ) THEN
     INFOG = -1
  ELSE IF ( N .LT. 0 ) THEN
     INFOG = -2
  ELSE IF ( LDA .LT. MAX( 1, N ) ) THEN
     INFOG = -4
  END IF

  NRANK = N
  IF ( INFO .GT. 0 ) NRANK = INFO - 1
  IF ( NRANK .EQ. 0 ) RETURN

  K = 1
10 CONTINUE
  IF ( K .GT. NRANK ) GOTO 20
  IF ( IPIV( K ) .GT. 0 ) THEN
     D = G( K, K )
     IF ( D .LT. ZERO ) THEN
        JVEC( K ) = -1
     ELSE
        JVEC( K ) = 1
     END IF
     D = SQRT( ABS( D ) )
     IF ( K .LT. N ) THEN
        CALL XSCAL( N-K, D, G( K+1, K ), 1 )
     END IF
     G( K, K ) = D
     K = K + 1
  ELSE IF ( IPIV( K ) .LT. 0 ) THEN
     AA = G( K, K )
     BB = G( K+1, K )
     CC = G( K+1, K+1 )
     CALL XLAEV2( AA, BB, CC, D, E, CS1, SN1 )
     IF ( D .LT. ZERO ) THEN
        JVEC( K ) = -1
        JVEC( K+1 ) = 1
     ELSE
        JVEC( K ) = 1
        JVEC( K+1 ) = -1
     END IF
     D = SQRT( ABS( D ) )
     E = SQRT( ABS( E ) )
     IF ( K+1 .LT. N ) THEN
        CALL XROT( N-K-1, G( K+2, K ), 1, G( K+2, K+1 ), 1, CS1, SN1 )
        CALL XSCAL( N-K-1, D, G( K+2, K ), 1 )
        CALL XSCAL( N-K-1, E, G( K+2, K+1 ), 1 )
     END IF
     G( K, K ) = D * CS1
     G( K+1, K+1 ) = E * CS1
     G( K, K+1 ) = -E * SN1
     G( K+1, K ) = D * SN1
     K = K + 2
  END IF
  GOTO 10

20 CONTINUE
  !     Make P*H*P' = U*J*U'   or  P*H*P' = L*J*L'

  !     Create the factor G such that  P*A*P' = G*J*G', where G is
  !     lower block triangular matrix, and J is as above.

  !     Empty the the the array G above the second superdiagonal
  NR2 = NRANK - 2
  IF ( NR2 .GT. 0 ) CALL XLASET( 'U', NR2, NR2, ZERO, ZERO, G( 1, 3 ), LDA )

  !     Set the element G( 1, 2 ) or G( 2, 3 ) to zero if necessary
  IF ( IPIV( 1 ) .GT. 0 ) THEN
     IF ( NRANK .GT. 1 ) G( 1, 2 ) = ZERO
     !     K is the loop index, increasing from 2 or 3 to NRANK in steps
     !     of 1 or 2, depending on the size of the diagonal blocks.
     K = 2
  ELSE
     IF ( NRANK .GT. 2 ) G( 2, 3 ) = ZERO
     K = 3
  END IF

30 CONTINUE
  !     If K > NRANK, exit from loop.
  IF ( K .GT. NRANK ) GOTO 40

  IF ( IPIV( K ) .GT. 0 ) THEN
     !     1 x 1 diagonal block

     !     Interchange rows K and IPIV(K).
     KP = IPIV( K )
     IF ( KP .NE. K ) CALL XSWAP( K-1, G( K, 1 ), LDA, G( KP, 1 ), LDA )

     !     Set the element G( k, k+1 ) of the upper triangle to zero
     IF ( NRANK .GT. K ) G( K, K+1 ) = ZERO
     K = K + 1
  ELSE
     !     2 x 2 diagonal block

     !     NO Interchange for row K.
     !     Interchange rows K+1 and -IPIV(K+1).
     KP = -IPIV( K+1 )
     IF ( KP .NE. K+1 ) CALL XSWAP( K-1, G( K+1, 1 ), LDA, G( KP, 1 ), LDA )

     !     Set the element G( k+1, k+2 ) of the upper triangle to zero
     IF ( NRANK .GT. K+1 ) G( K+1, K+2 ) = ZERO
     K = K + 2
  END IF
  GOTO 30

40 CONTINUE

END SUBROUTINE XGJGT
