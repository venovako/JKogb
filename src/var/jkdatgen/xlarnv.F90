SUBROUTINE XLARNV(IDIST, ISEED, N, X)
  IMPLICIT NONE
#include "wp.F"

!  -- LAPACK auxiliary routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

!     .. Scalar Arguments ..
  INTEGER, INTENT(IN) :: IDIST, N

!     .. Array Arguments ..
  INTEGER, INTENT(INOUT) :: ISEED(4)
  REAL(KIND=WP), INTENT(OUT) :: X(*)

!  Purpose
!  =======
!
!  XLARNV returns a vector of n random real numbers from a uniform or
!  normal distribution.
!
!  Arguments
!  =========
!
!  IDIST   (input) INTEGER
!          Specifies the distribution of the random numbers:
!          = 1:  uniform (0,1)
!          = 2:  uniform (-1,1)
!          = 3:  normal (0,1)
!
!  ISEED   (input/output) INTEGER array, dimension (4)
!          On entry, the seed of the random number generator; the array
!          elements must be between 0 and 4095, and ISEED(4) must be
!          odd.
!          On exit, the seed is updated.
!
!  N       (input) INTEGER
!          The number of random numbers to be generated.
!
!  X       (output) REAL(KIND=WP) array, dimension (N)
!          The generated random numbers.
!
!  Further Details
!  ===============
!
!  This routine calls the auxiliary routine XLARUV to generate random
!  real numbers from a uniform (0,1) distribution, in batches of up to
!  128 using vectorisable code. The Box-Muller method is used to
!  transform numbers from a uniform to a normal distribution.
!
!  =====================================================================

!     .. Parameters ..
  REAL(KIND=WP), PARAMETER :: ONE = +1.0E+0_WP, TWO = +2.0E+0_WP
  INTEGER, PARAMETER :: LV = 128
  REAL(KIND=WP), PARAMETER :: TWOPI = +6.2831853071795864769252867663E+0_WP

!     .. Local Scalars ..
  INTEGER :: I, IL, IL2, IV

!     .. Local Arrays ..
  REAL(KIND=WP) :: U(LV)

!     .. External Subroutines ..
  EXTERNAL :: XLARUV

!     .. Executable Statements ..

  DO IV = 1, N, LV / 2
     IL = MIN(LV / 2, N-IV+1)
     IF (IDIST .EQ. 3) THEN
        IL2 = 2*IL
     ELSE
        IL2 = IL
     END IF

!        Call XLARUV to generate IL2 numbers from a uniform (0,1)
!        distribution (IL2 <= LV)

     CALL XLARUV(ISEED, IL2, U)

     IF (IDIST .EQ. 1) THEN

!           Copy generated numbers

        DO I = 1, IL
           X( IV+I-1 ) = U( I )
        END DO
     ELSE IF (IDIST .EQ. 2) THEN

!           Convert generated numbers to uniform (-1,1) distribution

        DO I = 1, IL
           X( IV+I-1 ) = TWO*U( I ) - ONE
        END DO
     ELSE IF (IDIST .EQ. 3) THEN

!           Convert generated numbers to normal (0,1) distribution

        DO I = 1, IL
           X( IV+I-1 ) = SQRT(-TWO*LOG(U( 2*I-1 )))  * COS(TWOPI*U( 2*I ))
        END DO
     END IF

  END DO

END SUBROUTINE XLARNV
