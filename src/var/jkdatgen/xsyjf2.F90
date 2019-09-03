SUBROUTINE XSYJF2( UPLO, N, A, LDA, DIAGJ, IPIV, INFO )

  !     Modified by Singers, January 8, 2006.
  !
  !     Proposal by Ivan Slapnicar
  !     University of Split, Croatia
  !     slap@split.fesb.hr
  !     December 15, 1993
  !
  !     Purpose
  !     =======
  !
  !     XSYJF2 computes the factorization of a real symmetric matrix A using
  !     a modification of the Bunch-Kaufman diagonal pivoting method:
  !
  !     A = U*J*U'  or  A = L*J*L'
  !
  !     where U (or L) is a product of permutation and  upper (lower) block
  !     triangular matrices with 1-by-1 and 2-by-2 diagonal blocks, U' is
  !     the transpose of U, and J is diagonal with diagonal elements equal
  !     to 1 or -1 (see below for further details).
  !
  !     This is the unblocked version of the algorithm, calling Level 2 BLAS.
  !
  !     Arguments
  !     =========
  !
  !     UPLO    (input) CHARACTER*1
  !     Specifies whether the upper or lower triangular part of the
  !     symmetric matrix A is stored:
  !     = 'U':  Upper triangular
  !     = 'L':  Lower triangular
  !
  !     N       (input) INTEGER
  !     The order of the matrix A.  N >= 0.
  !
  !     A       (input/output) REAL(KIND=WP) array, dimension (LDA,N)
  !     On entry, the symmetric matrix A.  If UPLO = 'U', the leading
  !     N-by-N upper triangular part of A contains the upper
  !     triangular part of the matrix A. If UPLO = 'L', the
  !     leading N-by-N lower triangular part of A contains the lower
  !     triangular part of the matrix A.
  !
  !     On exit, the multipliers used to obtain the block triangular
  !     factor U or L (see below for further details).
  !
  !     LDA     (input) INTEGER
  !     The leading dimension of the array A.  LDA >= max(1,N).
  !
  !     DIAGJ   (output) REAL(KIND=WP) array, dimension(N)
  !     The diagonal of the matrix J.
  !
  !     IPIV    (output) INTEGER array, dimension (N)
  !     Details of the interchanges and the block structure of D.
  !     If IPIV(k) > 0, then rows and columns k and IPIV(k) were
  !     interchanged and in the k-th step a 1-by-1 diagonal block
  !     was used (see below for further details).
  !     If UPLO = 'U' and IPIV(k) < 0 and IPIV(k-1) < 0, then rows
  !     and columns k and -IPIV(k) and k-1 and -IPIV(k-1) were
  !     interchanged, and in the k-th step a 2-by-2 diagonal block
  !     was used.
  !     If UPLO = 'L' and IPIV(k) < 0 and IPIV(k+1) < 0, then rows
  !     and columns k and -IPIV(k) and k+1 and -IPIV(k+1) were
  !     interchanged, and in the k-th step a 2-by-2 diagonal block
  !     was used.
  !
  !     INFO    (output) INTEGER
  !     = 0: successful exit
  !     < 0: if INFO = -k, the k-th argument had an illegal value
  !     > 0: if UPLO = 'U' and INFO = k, then in the k-th step the
  !     leading submatrix A(1:k,1:k) was exactly zero, and the
  !     rank of A equals n-k;
  !     if UPLO = 'L' and INFO = k, then in the k-th step the
  !     trailing submatrix A(k:n,k:n) was exactly zero, and the
  !     rank of A equals k-1;
  !
  !     Further Details
  !     ===============
  !
  !     If UPLO = 'U', then A = U*J*U', where
  !     U = P(n)*U(n)* ... *P(k)U(k)* ...,
  !     i.e., U is a product of terms P(k)*U(k), where k decreases from n to
  !     1 in steps of 1 or 2, and J is a diagonal matrix with diagonal
  !     elements equal to 1 or -1. P(k) is a permutation matrix as defined
  !     by IPIV(k), and U(k) is a upper block triangular matrix, such that
  !
  !             (   I    v    0   )   k-s
  !     U(k) =  (   0    d    0   )   s
  !             (   0    0    I   )   n-k
  !                k-s   s   n-k
  !
  !     Here s = 1 or 2, and ( v ) overwrites A(1:k,k-s+1:k).
  !                          ( d )
  !     
  !     If UPLO = 'L', then A = L*J*L', where
  !     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
  !     i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
  !     n in steps of 1 or 2, and J is a diagonal matrix with diagonal
  !     elements equal to 1 or -1. P(k) is a permutation matrix as defined
  !     by IPIV(k), and L(k) is a lower block triangular matrix, such that
  !
  !             (   I    0     0   )  k-1
  !     L(k) =  (   0    d     0   )  s
  !             (   0    v     I   )  n-k-s+1
  !                k-1   s  n-k-s+1
  !
  !     Here s = 1 or 2, and ( d ) overwrites A(k:n,k:k+s-1).
  !                          ( v )
  !
  !
  !     The following comment is related to a comment in the eigenreduction
  !     routine SSYEJA: the output format
  !
  !     U = P(n)*U(n)* ... *P(k)U(k)* ...,  or
  !     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
  !
  !     is taken from the Lapack routine SSYTF2. For our purpose it might be
  !     more convenient to have the output in the form
  !     P*H*P' = U*J*U'   or  P*H*P' = L*J*L'
  !     as in the original Bunch and Parlett paper from Wilkinson, Reinsch book,
  !     and in the former routine SGJGT. Then the first backward permutation
  !     in SSYEJA would not be necessary. Thus, the final choice is an open
  !     question.

  IMPLICIT NONE
#include "wp.F"

  CHARACTER(LEN=1), INTENT(IN) :: UPLO
  INTEGER, INTENT(IN) :: N
  REAL(KIND=WP), DIMENSION(LDA, N), INTENT(INOUT) :: A
  INTEGER, INTENT(IN) :: LDA
  REAL(KIND=WP), DIMENSION(N), INTENT(OUT) :: DIAGJ
  INTEGER, DIMENSION(N), INTENT(OUT) :: IPIV
  INTEGER, INTENT(OUT) :: INFO

  REAL(KIND=WP), PARAMETER ::   ZERO =  +0.0E+0_WP
  REAL(KIND=WP), PARAMETER ::    ONE =  +1.0E+0_WP
  REAL(KIND=WP), PARAMETER ::  EIGHT =  +8.0E+0_WP
  REAL(KIND=WP), PARAMETER :: SEVTEN = +17.0E+0_WP

  LOGICAL :: UPPER
  INTEGER :: I, J, IMAX, JMAX, K, KDIAG, KP, KSTEP
  REAL(KIND=WP) :: ALPHA, C, DIAMAX, OFFMAX, R1, R2, S, T, TEMP

  INTEGER, EXTERNAL :: IXAMAX
  LOGICAL, EXTERNAL :: LSAME

  EXTERNAL :: XLAEV2, XROT, XSCAL, XSWAP, XSYR

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

  !     Initialize ALPHA for use in choosing pivot block size.
  ALPHA = ( ONE + SQRT( SEVTEN ) ) / EIGHT
  IMAX = 0
  JMAX = 0

  IF ( UPPER ) THEN
     !     Factorize A as U*J*U' using the upper triangle of A
     !
     !     K is the main loop index, decreasing from N to 1 in steps of 1 or 2
     K = N
10   CONTINUE
     !     If K < 1, exit from loop
     IF ( K .LT. 1 ) GOTO 70

     !     Determine rows and columns to be interchanged and whether
     !     a 1-by-1 or 2-by-2 pivot block will be used

     !     KDIAG is the index of the largest diagonal element, and
     !     DIAMAX is its absolute value
     KDIAG = IXAMAX( K, A( 1, 1 ), LDA + 1 )
     DIAMAX = ABS( A( KDIAG, KDIAG ) )

     !     IMAX and JMAX are the row- and column-indices of the largest
     !     off-diagonal element, and OFFMAX is its absolute value
     OFFMAX = ZERO
     IF ( K .GT. 1 ) THEN
        DO J = 2, K
           DO I = 1, J - 1
              TEMP = ABS( A( I, J ) )
              IF ( TEMP .GT. OFFMAX ) THEN
                 OFFMAX = TEMP
                 IMAX = I
                 JMAX = J
              END IF
           END DO
        END DO
     END IF

     IF ( MAX( DIAMAX, OFFMAX ) .EQ. ZERO ) THEN
        !     The rest of the matrix is zero: set INFO and return
        INFO = K
        GOTO 70
     END IF

     IF ( DIAMAX .GE. ALPHA * OFFMAX ) THEN
        !     Use 1-by-1 pivot block
        KSTEP = 1
        KP = KDIAG
     ELSE
        !     Use 2-by-2 pivot block
        KSTEP = 2
        KP = JMAX
     END IF

     IF ( KP .NE. K ) THEN
        !     Interchange rows and columns K and KP in the leading
        !     submatrix A(1:k,1:k)
        CALL XSWAP( KP-1, A( 1, K ), 1, A( 1, KP ), 1 )
        CALL XSWAP( K-KP-1, A( KP+1, K ), 1, A( KP, KP+1 ), LDA )
        T = A( K, K )
        A( K, K ) = A( KP, KP )
        A( KP, KP ) = T
     END IF

     IF ( ( KSTEP .EQ. 2 ) .AND. ( IMAX .NE. K-1 ) ) THEN
        !     Interchange rows and columns K-1 and IMAX in the leading
        !     submatrix A(1:k,1:k)
        CALL XSWAP( IMAX-1, A( 1, K-1 ), 1, A( 1, IMAX ), 1 )
        CALL XSWAP( K-IMAX-2, A( IMAX+1, K-1 ), 1, A( IMAX, IMAX+1 ), LDA )
        T = A( K-1, K-1 )
        A( K-1, K-1 ) = A( IMAX, IMAX )
        A( IMAX, IMAX ) = T
        T = A( IMAX, K )
        A( IMAX, K ) = A( K-1, K )
        A( K-1, K ) = T
     END IF

     !     Update the leading submatrix
     IF ( KSTEP .EQ. 1 ) THEN
        !     1-by-1 pivot block D(k): column k now holds
        !
        !     W(1:k-1,k) = U(1:k-1,k) * sqrt( abs(D(k)) ),
        !     W(k,k) = U(k,k) * sqrt( abs(D(k)) ) * J(k),
        !
        !     where U(1:k,k) is the k-th column of U, and
        !     J(k) = sign( D(k) )
        !
        !     Perform a rank-1 update of A(1:k-1,1:k-1) as
        !
        !     A := A - U(k)*J(k,k)*U(k)' = A - W(k)*1/D(k)*W(k)'
        R1 = ONE / A( K, K )
        IF ( K .GT. 1 ) CALL XSYR( UPLO, K-1, -R1, A( 1, K ), 1, A, LDA )

        !     Compute the k-th diagonal element of the matrix J
        IF ( R1 .GE. ZERO ) THEN
           DIAGJ( K ) = ONE
        ELSE
           DIAGJ( K ) = -ONE
        END IF

        !     Store U(k) in column k
        R1 = SQRT( ABS( R1 ) )
        CALL XSCAL( K-1, R1, A( 1, K ), 1 )
        A( K, K ) = DIAGJ( K ) / R1
     ELSE
        !     2-by-2 pivot block D(k): let
        !
        !     D(k) = Q(k)**T * X(k) * Q(k)
        !
        !     be the eigendecomposition of D(k), X(k) = diag(R1,R2).
        !     Columns k and k-1 now hold
        !
        !     ( W(1:k-2,k-1) W(1:k-2,k) ) =
        !
        !     ( U(1:k-2,k-1) U(1:k-2,k) )*sqrt(abs(X(k)))*Q(k)**T,
        !
        !     W(k-1:k,k-1:k) =
        !
        !     U(k-1:k,k-1:k)*Q(k)*inv(sqrt(abs(X(k))))*J(k),
        !
        !     where U(k) and U(k-1) are the k-th and (k-1)-th columns
        !     of U, and J(k) = diag( sign(R1), sign(R2) ).
        !
        !     Perform a rank-2 update of A(1:k-2,1:k-2) as
        !
        !     A := A - ( U(k-1) U(k) )*J(k)*( U(k-1) U(k) )'
        !     = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )'
        !
        !     Convert this to two rank-1 updates by using the eigen-
        !     decomposition of D(k)
        CALL XLAEV2( A( K-1, K-1 ), A( K-1, K ), A( K, K ), R1, R2, C, S )
        R1 = ONE / R1
        R2 = ONE / R2
        CALL XROT( K-2, A( 1, K-1 ), 1, A( 1, K ), 1, C, S )
        CALL XSYR( UPLO, K-2, -R1, A( 1, K-1 ), 1, A, LDA )
        CALL XSYR( UPLO, K-2, -R2, A( 1, K ), 1, A, LDA )

        !     Compute the (k-1)-st and k-th diagonal element of the matrix J
        IF ( R1 .GE. ZERO ) THEN
           DIAGJ( K-1 ) = ONE
           DIAGJ( K ) = -ONE
        ELSE
           DIAGJ( K-1 ) = -ONE
           DIAGJ( K ) = ONE
        END IF

        !     Store U(k) and U(k-1) in columns k and k-1
        R1 = SQRT( ABS( R1 ) )
        R2 = SQRT( ABS( R2 ) )
        CALL XSCAL( K-2, R1, A( 1, K-1 ), 1 )
        CALL XSCAL( K-2, R2, A( 1, K ), 1 )

        R1 = DIAGJ( K-1 ) / R1
        R2 = DIAGJ( K ) / R2
        A( K-1, K-1 ) = C * R1
        A( K, K-1 ) = S * R1
        A( K-1, K ) = -S * R2
        A( K, K ) = C * R2
     END IF

     !     Store details of the interchanges in IPIV
     IF ( KSTEP .EQ. 1 ) THEN
        IPIV( K ) = KP
     ELSE
        IPIV( K ) = -KP
        IPIV( K-1 ) = -IMAX
     END IF

     !     Decrease K and return to the start of the main loop
     K = K - KSTEP
     GOTO 10
  ELSE
     !     Factorize A as L*D*L' using the lower triangle of A
     !
     !     K is the main loop index, increasing from 1 to N in steps of 1 or 2
     K = 1
40   CONTINUE
     !     If K > N, exit from loop
     IF ( K .GT. N ) GOTO 70

     !     Determine rows and columns to be interchanged and whether
     !     a 1-by-1 or 2-by-2 pivot block will be used

     !     KDIAG is the index of the largest diagonal element, and
     !     DIAMAX is its absolute value
     KDIAG = IXAMAX( N-K+1, A( K, K ), LDA + 1 ) + K - 1
     DIAMAX = ABS( A( KDIAG, KDIAG ) )

     !     IMAX and JMAX are the row- and column-indices of the largest
     !     off-diagonal element, and OFFMAX is its absolute value
     OFFMAX = ZERO
     IF ( K .LT. N ) THEN
        DO J = K, N - 1
           DO I = J + 1, N
              TEMP = ABS( A( I, J ) )
              IF ( TEMP .GT. OFFMAX ) THEN
                 OFFMAX = TEMP
                 IMAX = I
                 JMAX = J
              END IF
           END DO
        END DO
     END IF

     IF ( MAX( DIAMAX, OFFMAX ) .EQ. ZERO ) THEN
        !     The rest of the matrix is zero: set INFO and return
        INFO = K
        GOTO 70
     END IF

     IF ( DIAMAX .GE. ALPHA * OFFMAX ) THEN
        !     Use 1-by-1 pivot block
        KSTEP = 1
        KP = KDIAG
     ELSE
        !     Use 2-by-2 pivot block
        KSTEP = 2
        KP = JMAX
     END IF

     IF ( KP .NE. K ) THEN
        !     Interchange rows and columns K and KP in the trailing
        !     submatrix A(k:n,k:n)
        IF ( KP .LT. N ) CALL XSWAP( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 )
        CALL XSWAP( KP-K-1, A( K+1, K ), 1, A( KP, K+1 ), LDA )
        T = A( K, K )
        A( K, K ) = A( KP, KP )
        A( KP, KP ) = T
     END IF

     IF ( ( KSTEP .EQ. 2 ) .AND. ( IMAX .NE. K+1 ) ) THEN
        !     Interchange rows and columns K+1 and IMAX in the trailing
        !     submatrix A(k:n,k:n)
        IF ( IMAX .LT. N ) CALL XSWAP( N-IMAX, A( IMAX+1, K+1 ), 1, A( IMAX+1, IMAX ), 1 )
        CALL XSWAP( IMAX-K-2, A( K+2, K+1 ), 1, A( IMAX, K+2 ), LDA )
        T = A( K+1, K+1 )
        A( K+1, K+1 ) = A( IMAX, IMAX )
        A( IMAX, IMAX ) = T
        T = A( IMAX, K )
        A( IMAX, K ) = A( K+1, K )
        A( K+1, K ) = T
     END IF

     !     Update the trailing submatrix
     IF ( KSTEP .EQ. 1 ) THEN
        !     1-by-1 pivot block D(k): column k now holds
        !
        !     W(k+1:n,k) = U(k+1:n,k) * sqrt( abs(D(k)) ),
        !     W(k,k) = U(k,k) * sqrt( abs(D(k)) ) * J(k),
        !
        !     where U(k:n,k) is the k-th column of U, and
        !     J(k) = sign( D(k) )
        !
        !     Perform a rank-1 update of A(k+1:n,k+1:n) as
        !
        !     A := A - U(k)*J(k,k)*U(k)' = A - W(k)*1/D(k)*W(k)'
        R1 = ONE / A( K, K )
        IF ( K .LT. N ) CALL XSYR( UPLO, N-K, -R1, A( K+1, K ), 1, A( K+1, K+1 ), LDA )

        !     Compute the k-th diagonal element of the matrix J
        IF ( R1 .GE. ZERO ) THEN
           DIAGJ( K ) = ONE
        ELSE
           DIAGJ( K ) = -ONE
        END IF

        !     Store U(k) in column k
        R1 = SQRT( ABS( R1 ) )
        !     Singer modif: IF k < n test added!
        IF ( K .LT. N ) CALL XSCAL( N-K, R1, A( K+1, K ), 1 )
        A( K, K ) = DIAGJ( K ) / R1
     ELSE
        !     2-by-2 pivot block D(k): let
        !
        !     D(k) = Q(k)**T * X(k) * Q(k)
        !
        !     be the eigendecomposition of D(k), X(k) = diag(R1,R2).
        !     Columns k and k-1 now hold
        !
        !     ( W(k+2:n,k) W(k+2:n,k+1) ) =
        !
        !     ( U(k+2:n,k) U(k+2:n,k+1) )*sqrt(abs(X(k)))*Q(k)**T,
        !
        !     W(k:k+1,k:k+1) =
        !
        !     U(k:k+1,k:k+1)*Q(k)*inv(sqrt(abs(X(k))))*J(k),
        !
        !     where U(k) and U(k+1) are the k-th and (k+1)-st columns
        !     of U, and J(k) = diag( sign(R1), sign(R2) ).
        !
        !     Perform a rank-2 update of A(k+2:n,k+2:n) as
        !
        !     A := A - ( U(k) U(k+1) )*J(k)*( U(k) U(k+1) )'
        !     = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )'
        !
        !     Convert this to two rank-1 updates by using the eigen-
        !     decomposition of D(k)
        CALL XLAEV2( A( K, K ), A( K+1, K ), A( K+1, K+1 ), R1, R2, C, S )
        R1 = ONE / R1
        R2 = ONE / R2
        CALL XROT( N-K-1, A( K+2, K ), 1, A( K+2, K+1 ), 1, C, S )
        CALL XSYR( UPLO, N-K-1, -R1, A( K+2, K ), 1, A( K+2, K+2 ), LDA )
        CALL XSYR( UPLO, N-K-1, -R2, A( K+2, K+1 ), 1, A( K+2, K+2 ), LDA )

        !     Compute the k-th and (k+1)-st diagonal element of the matrix J
        IF ( R1 .GE. ZERO ) THEN
           DIAGJ( K ) = ONE
           DIAGJ( K+1 ) = -ONE
        ELSE
           DIAGJ( K ) = -ONE
           DIAGJ( K+1 ) = ONE
        END IF

        !     Store U(k) and U(k+1) in columns k and k+1
        R1 = SQRT( ABS( R1 ) )
        R2 = SQRT( ABS( R2 ) )
        CALL XSCAL( N-K-1, R1, A( K+2, K ), 1 )
        CALL XSCAL( N-K-1, R2, A( K+2, K+1 ), 1 )

        R1 = DIAGJ( K ) / R1
        R2 = DIAGJ( K+1 ) / R2
        A( K, K ) = C * R1
        A( K, K+1 ) = -S * R2
        A( K+1, K ) = S * R1
        A( K+1, K+1 ) = C * R2
     END IF

     !     Store details of the interchanges in IPIV
     IF ( KSTEP .EQ. 1 ) THEN
        IPIV( K ) = KP
     ELSE
        IPIV( K ) = -KP
        IPIV( K+1 ) = -IMAX
     END IF

     !     Increase K and return to the start of the main loop
     K = K + KSTEP
     GOTO 40
  END IF

70 CONTINUE

END SUBROUTINE XSYJF2
