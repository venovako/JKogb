FUNCTION XLARND(IDIST, ISEED)
  IMPLICIT NONE
#include "wp.F"

!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..

!     .. Scalar Arguments ..
  INTEGER, INTENT(IN) :: IDIST

!     .. Array Arguments ..
  INTEGER, INTENT(INOUT) :: ISEED(4)

  REAL(KIND=WP) :: XLARND

!  Purpose
!  =======
!
!  XLARND returns a random real number from a uniform or normal
!  distribution.
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
!  Further Details
!  ===============
!
!  This routine calls the auxiliary routine XLARAN to generate a random
!  real number from a uniform (0,1) distribution. The Box-Muller method
!  is used to transform numbers from a uniform to a normal distribution.
!
!  =====================================================================

!     .. Parameters ..
  REAL(KIND=WP), PARAMETER :: ZERO = +0.0E+0_WP, ONE = +1.0E+0_WP, TWO = +2.0E+0_WP
  REAL(KIND=WP), PARAMETER :: TWOPI = +6.2831853071795864769252867663E+0_WP

!     .. Local Scalars ..
  REAL(KIND=WP) :: T1, T2

!     .. External Functions ..
  REAL(KIND=WP), EXTERNAL :: XLARAN

!     .. Executable Statements ..

!     Generate a real random number from a uniform (0,1) distribution

  T1 = XLARAN(ISEED)

  IF (IDIST .EQ. 1) THEN
!        uniform (0,1)
     XLARND = T1
  ELSE IF (IDIST .EQ. 2) THEN
!        uniform (-1,1)
     XLARND = TWO*T1 - ONE
  ELSE IF (IDIST .EQ. 3) THEN
!        normal (0,1)
     T2 = XLARAN(ISEED)
     XLARND = SQRT(-TWO*LOG(T1)) * COS(TWOPI*T2)
  ELSE
!        error
     XLARND = ZERO
  END IF

END FUNCTION XLARND
