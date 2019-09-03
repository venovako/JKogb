SUBROUTINE XSYBPC( UPLO, N, A, LDA, NRANK, DIAGJ, IPIV, JVEC, INFO )

  !     Modified by Singers, January 8, 2006.
  !
  !     Purpose
  !     =======
  !
  !     XSYBPC computes the modified Bunch-Parlett factorization of a N-by-N
  !     real symmetric matrix A
  !
  !     A = G * J * transpose( G ).
  !
  !     Arguments
  !     =========
  !
  !     UPLO    (input) CHARACTER*1
  !     = 'U':  Upper triangle of A is stored;
  !     = 'L':  Lower triangle of A is stored.
  !
  !     N       (input) INTEGER
  !     Order of the input matrix A, N >= 0.
  !
  !     A       (input/output) REAL(KIND=WP) array, dimension (LDA,N)
  !     On entry, the symmetric matrix A.  If UPLO = 'U', the leading
  !     N-by-N upper triangular part of A contains the upper
  !     triangular part of the matrix A. If UPLO = 'L', the
  !     leading N-by-N lower triangular part of A contains the lower
  !     triangular part of the matrix A.
  !
  !     On exit, A contains the matrix G.
  !
  !     LDA     (input) INTEGER
  !     The leading dimension of the array A.  LDA >= max(1,N).
  !
  !     NRANK   (output) INTEGER
  !     Contains the rank of A.
  !
  !     DIAGJ   (workspace) REAL(KIND=WP) array, dimension(N)
  !     Contains the signs of the eigenvalues.
  !
  !     IPIV    (workspace) INTEGER array, dimension(N)
  !
  !     JVEC    (output) INTEGER array, dimension(N)
  !     Contains the diagonal of the matrix J; JVEC( I ) = 1 or -1.
  !     If NRANK < N, only the first NRANK values are set.
  !
  !     INFO    (output) INTEGER
  !     = 0:  successful exit - the first N columns of the array A
  !     contain the matrix G (full column rank).
  !     < 0:  if INFO = -i, with 0 < i <= 1000, the i-th argument
  !     had an illegal value;
  !     > 0:  some eigenvalues are zero and INFO specifies the
  !     rank of A.

  IMPLICIT NONE
#include "wp.F"

  CHARACTER(LEN=1), INTENT(IN) :: UPLO
  INTEGER, INTENT(IN) :: N
  REAL(KIND=WP), DIMENSION(LDA, N), INTENT(INOUT) :: A
  INTEGER, INTENT(IN) :: LDA
  INTEGER, INTENT(OUT) :: NRANK
  REAL(KIND=WP), DIMENSION(N), INTENT(OUT) :: DIAGJ
  INTEGER, DIMENSION(N), INTENT(OUT) :: IPIV
  INTEGER, DIMENSION(N), INTENT(OUT) :: JVEC
  INTEGER, INTENT(OUT) :: INFO

  REAL(KIND=WP), PARAMETER :: ZERO = +0.0E+0_WP

  LOGICAL :: UPPER
  INTEGER :: I, INFOD
  INTEGER :: K, KP, NR2

  LOGICAL, EXTERNAL :: LSAME
  EXTERNAL :: XCOPY, XLASET, XSWAP, XSYJF2

  !     Test the input arguments
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )

  IF ( .NOT. UPPER .AND. .NOT. LSAME( UPLO, 'L' ) ) THEN
     INFO = -1
  ELSE IF ( N .LT. 0 ) THEN
     INFO = -2
  ELSE IF ( LDA .LT. MAX( 1, N ) ) THEN
     INFO = -4
  END IF

  IF ( INFO .NE. 0 ) RETURN

  !     Quick return, if possible.
  IF ( N .EQ. 0 ) THEN
     NRANK = 0
     RETURN
  END IF

  !     Compute the factorization from the lower triangle of A.
  IF ( UPPER ) THEN
     !     Copy the upper triangle to the lower triangle
     DO I = 1, N - 1
        CALL XCOPY( N-I, A( I, I+1 ), LDA, A( I+1, I ), 1 )
     END DO
  END IF

  !     Compute the factorization A = L*J*L', where L is a product
  !     of permutation and lower block triangular matrices with
  !     1-by-1 and 2-by-2 diagonal blocks, L' is the transpose of L,
  !     and J is diagonal with diagonal elements equal to 1 or -1.
  CALL XSYJF2( 'L', N, A, LDA, DIAGJ, IPIV, INFOD )

  !     Set NRANK to the rank of A
  IF ( INFOD .EQ. 0 ) THEN
     NRANK = N
  ELSE
     NRANK = INFOD - 1
  END IF

  !     Create the factor G such that  P*A*P' = G*J*G', where G is
  !     lower block triangular matrix, and J is as above.

  !     Empty the the the array A above the second superdiagonal
  NR2 = NRANK - 2
  IF ( NR2 .GT. 0 ) CALL XLASET( 'U', NR2, NR2, ZERO, ZERO, A( 1, 3 ), LDA )

  !     Set the element A( 1, 2 ) or A( 2, 3 ) to zero if necessary
  IF ( IPIV( 1 ) .GT. 0 ) THEN
     IF ( NRANK .GT. 1 ) A( 1, 2 ) = ZERO
     !     K is the loop index, increasing from 2 or 3 to NRANK in steps
     !     of 1 or 2, depending on the size of the diagonal blocks.
     K = 2
  ELSE
     IF ( NRANK .GT. 2 ) A( 2, 3 ) = ZERO
     K = 3
  END IF

20 CONTINUE
  !     If K > NRANK, exit from loop.
  IF ( K .GT. NRANK ) GOTO 30

  IF ( IPIV( K ) .GT. 0 ) THEN
     !     1 x 1 diagonal block
     !
     !     Interchange rows K and IPIV(K).
     KP = IPIV( K )
     IF ( KP .NE. K ) CALL XSWAP( K-1, A( K, 1 ), LDA, A( KP, 1 ), LDA )

     !     Set the element A( k, k+1 ) of the upper triangle to zero
     IF ( NRANK .GT. K ) A( K, K+1 ) = ZERO
     K = K + 1
  ELSE
     !     2 x 2 diagonal block
     !
     !     Interchange rows K and -IPIV(K).
     KP = -IPIV( K )
     IF ( KP .NE. K ) CALL XSWAP( K-1, A( K, 1 ), LDA, A( KP, 1 ), LDA )

     !     Interchange rows K+1 and -IPIV(K+1).
     KP = -IPIV( K+1 )
     IF ( KP .NE. K+1 ) CALL XSWAP( K-1, A( K+1, 1 ), LDA, A( KP, 1 ), LDA )

     !     Set the element A( K+1, K+2 ) of the upper triangle to zero
     IF ( NRANK .GT. K+1 ) A( K+1, K+2 ) = ZERO
     K = K + 2
  END IF
  GOTO 20

30 CONTINUE
  DO I = 1, NRANK
     JVEC( I ) = INT( DIAGJ( I ) )
  END DO

  !     If column rank defect occured, set INFO = RANK
  IF ( NRANK .LT. N ) INFO = NRANK

END SUBROUTINE XSYBPC
