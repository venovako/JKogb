real(kind=16) function xlarnd(idist, iseed)

  implicit none

!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..

!     .. Scalar Arguments ..
  integer, intent(in) :: idist

!     .. Array Arguments ..
  integer, intent(inout) :: iseed(4)

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
  real(kind=16), parameter :: zero = +0.0e+0_16, one = +1.0e+0_16, two = +2.0e+0_16
  real(kind=16), parameter :: twopi = +6.2831853071795864769252867663e+0_16

!     .. Local Scalars ..
  real(kind=16) :: t1, t2

!     .. External Functions ..
  real(kind=16), external :: xlaran

!     .. Executable Statements ..

!     Generate a real random number from a uniform (0,1) distribution

  t1 = xlaran(iseed)

  if (idist .eq. 1) then
!        uniform (0,1)
     xlarnd = t1
  else if (idist .eq. 2) then
!        uniform (-1,1)
     xlarnd = two*t1 - one
  else if (idist .eq. 3) then
!        normal (0,1)
     t2 = xlaran(iseed)
     xlarnd = sqrt(-two*log(t1)) * cos(twopi*t2)
  else
!        error
     xlarnd = zero
  end if

end function xlarnd
