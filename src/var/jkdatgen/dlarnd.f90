real(kind=8) function dlarnd( idist, iseed )

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
!  DLARND returns a random real number from a uniform or normal
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
!  This routine calls the auxiliary routine DLARAN to generate a random
!  real number from a uniform (0,1) distribution. The Box-Muller method
!  is used to transform numbers from a uniform to a normal distribution.
!
!  =====================================================================

!     .. Parameters ..
  real(kind=8), parameter :: zero = +0.0e+0_8, one = +1.0e+0_8, two = +2.0e+0_8
  real(kind=8), parameter :: twopi = +6.2831853071795864769252867663e+0_8

!     .. Local Scalars ..
  real(kind=8) :: t1, t2

!     .. External Functions ..
  real(kind=8), external :: dlaran

!     .. Executable Statements ..

!     Generate a real random number from a uniform (0,1) distribution

  t1 = dlaran( iseed )

  if (idist .eq. 1) then

!        uniform (0,1)

     dlarnd = t1
  else if (idist .eq. 2) then

!        uniform (-1,1)

     dlarnd = two*t1 - one
  else if (idist .eq. 3) then

!        normal (0,1)

     t2 = dlaran( iseed )
     dlarnd = sqrt(-two*log(t1)) * cos(twopi*t2)
  else

!        error

     dlarnd = zero
  end if

end function dlarnd
