SUBROUTINE XLAGSY(N, K, D, A, LDA, ISEED, WORK, INFO)
  IMPLICIT NONE
#include "wp.F"

!  -- LAPACK auxiliary test routine (version 3.1)
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..

!     .. Scalar Arguments ..
  INTEGER, INTENT(OUT) :: INFO
  INTEGER, INTENT(IN) :: K, LDA, N

!     .. Array Arguments ..
  INTEGER, INTENT(INOUT) :: ISEED(4)
  REAL(KIND=WP), INTENT(OUT) :: A(LDA,*)
  REAL(KIND=WP), INTENT(IN) :: D(*)
  REAL(KIND=WP), INTENT(OUT) :: WORK(*)

!  Purpose
!  =======
!
!  XLAGSY generates a real symmetric matrix A, by pre- and post-
!  multiplying a real diagonal matrix D with a random orthogonal matrix:
!  A = U*D*U'. The semi-bandwidth may then be reduced to k by additional
!  orthogonal transformations.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  K       (input) INTEGER
!          The number of nonzero subdiagonals within the band of A.
!          0 <= K <= N-1.
!
!  D       (input) REAL(KIND=WP) array, dimension (N)
!          The diagonal elements of the diagonal matrix D.
!
!  A       (output) REAL(KIND=WP) array, dimension (LDA,N)
!          The generated n by n symmetric matrix A (the full matrix is
!          stored).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= N.
!
!  ISEED   (input/output) INTEGER array, dimension (4)
!          On entry, the seed of the random number generator; the array
!          elements must be between 0 and 4095, and ISEED(4) must be
!          odd.
!          On exit, the seed is updated.
!
!  WORK    (workspace) REAL(KIND=WP) array, dimension (2*N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================

!     .. Parameters ..
  REAL(KIND=WP), PARAMETER :: ZERO = +0.0E+0_WP, ONE = +1.0E+0_WP, HALF = +0.5E+0_WP

!     .. Local Scalars ..
  INTEGER :: I, J
  REAL(KIND=WP) :: ALPHA, TAU, WA, WB, WN

!     .. External Subroutines ..
  EXTERNAL :: XAXPY, XGEMV, XGER, XLARNV, XSCAL, XSYMV, XSYR2, XERBLA

!     .. External Functions ..
  REAL(KIND=WP), EXTERNAL :: XDOT, XNRM2

!     .. Executable Statements ..

!     Test the input arguments

  INFO = 0
  IF (N .LT. 0) THEN
     INFO = -1
  ELSE IF (K .LT. 0 .OR. K .GT. N-1) THEN
     INFO = -2
  ELSE IF (LDA .LT. MAX(1,N)) THEN
     INFO = -5
  END IF
  IF (INFO .LT. 0) THEN
     CALL XERBLA('XLAGSY', -INFO)
     RETURN
  END IF

!     initialize lower triangle of A to diagonal matrix

  DO J = 1, N
     DO I = J + 1, N
        A( I, J ) = ZERO
     END DO
  END DO
  DO I = 1, N
     A( I, I ) = D( I )
  END DO

!     Generate lower triangle of symmetric matrix

  DO I = N - 1, 1, -1

!        generate random reflection

     CALL XLARNV(3, ISEED, N-I+1, WORK)
     WN = XNRM2(N-I+1, WORK, 1)
     WA = SIGN(WN, WORK( 1 ))
     IF (WN .EQ. ZERO) THEN
        TAU = ZERO
     ELSE
        WB = WORK( 1 ) + WA
        CALL XSCAL(N-I, ONE / WB, WORK( 2 ), 1)
        WORK( 1 ) = ONE
        TAU = WB / WA
     END IF

!        apply random reflection to A(i:n,i:n) from the left
!        and the right

!        compute  y := tau * A * u

     CALL XSYMV('L', N-I+1, TAU, A( I, I ), LDA, WORK, 1, ZERO, WORK( N+1 ), 1)

!        compute  v := y - 1/2 * tau * ( y, u ) * u

     ALPHA = -HALF * TAU * XDOT( N-I+1, WORK( N+1 ), 1, WORK, 1 )
     CALL XAXPY(N-I+1, ALPHA, WORK, 1, WORK( N+1 ), 1)

!        apply the transformation as a rank-2 update to A(i:n,i:n)

     CALL XSYR2('L', N-I+1, -ONE, WORK, 1, WORK( N+1 ), 1, A( I, I ), LDA)
  END DO

!     Reduce number of subdiagonals to K

  DO I = 1, N - 1 - K

!        generate reflection to annihilate A(k+i+1:n,i)

     WN = XNRM2(N-K-I+1, A( K+I, I ), 1)
     WA = SIGN(WN, A( K+I, I ))
     IF (WN .EQ. ZERO) THEN
        TAU = ZERO
     ELSE
        WB = A( K+I, I ) + WA
        CALL XSCAL(N-K-I, ONE / WB, A( K+I+1, I ), 1)
        A( K+I, I ) = ONE
        TAU = WB / WA
     END IF

!        apply reflection to A(k+i:n,i+1:k+i-1) from the left

     CALL XGEMV('T', N-K-I+1, K-1, ONE, A( K+I, I+1 ), LDA, A( K+I, I ), 1, ZERO, WORK, 1)
     CALL XGER(N-K-I+1, K-1, -TAU, A( K+I, I ), 1, WORK, 1, A( K+I, I+1 ), LDA)

!        apply reflection to A(k+i:n,k+i:n) from the left and the right

!        compute  y := tau * A * u

     CALL XSYMV('L', N-K-I+1, TAU, A( K+I, K+I ), LDA, A( K+I, I ), 1, ZERO, WORK, 1)

!        compute  v := y - 1/2 * tau * ( y, u ) * u

     ALPHA = -HALF * TAU * XDOT(N-K-I+1, WORK, 1, A( K+I, I ), 1)
     CALL XAXPY(N-K-I+1, ALPHA, A( K+I, I ), 1, WORK, 1)

!        apply symmetric rank-2 update to A(k+i:n,k+i:n)

     CALL XSYR2('L', N-K-I+1, -ONE, A( K+I, I ), 1, WORK, 1, A( K+I, K+I ), LDA)

     A( K+I, I ) = -WA
     DO J = K + I + 1, N
        A( J, I ) = ZERO
     END DO
  END DO

!     Store full symmetric matrix

  DO J = 1, N
     DO I = J + 1, N
        A( J, I ) = A( I, J )
     END DO
  END DO

END SUBROUTINE XLAGSY
