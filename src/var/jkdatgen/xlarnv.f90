subroutine xlarnv(idist, iseed, n, x)

  implicit none

!  -- LAPACK auxiliary routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

!     .. Scalar Arguments ..
  integer, intent(in) :: idist, n

!     .. Array Arguments ..
  integer, intent(inout) :: iseed(4)
  real(kind=16), intent(out) :: x(*)

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
!  X       (output) REAL(KIND=16) array, dimension (N)
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
  real(kind=16), parameter :: one = +1.0e+0_16, two = +2.0e+0_16
  integer, parameter :: lv = 128
  real(kind=16), parameter :: twopi = +6.2831853071795864769252867663e+0_16

!     .. Local Scalars ..
  integer :: i, il, il2, iv

!     .. Local Arrays ..
  real(kind=16) :: u(lv)

!     .. External Subroutines ..
  external :: xlaruv

!     .. Executable Statements ..

  do iv = 1, n, lv / 2
     il = min(lv / 2, n-iv+1)
     if (idist .eq. 3) then
        il2 = 2*il
     else
        il2 = il
     end if

!        Call XLARUV to generate IL2 numbers from a uniform (0,1)
!        distribution (IL2 <= LV)

     call xlaruv(iseed, il2, u)

     if (idist .eq. 1) then

!           Copy generated numbers

        do i = 1, il
           x( iv+i-1 ) = u( i )
        end do
     else if (idist .eq. 2) then

!           Convert generated numbers to uniform (-1,1) distribution

        do i = 1, il
           x( iv+i-1 ) = two*u( i ) - one
        end do
     else if (idist .eq. 3) then

!           Convert generated numbers to normal (0,1) distribution

        do i = 1, il
           x( iv+i-1 ) = sqrt(-two*log(u( 2*i-1 )))  * cos(twopi*u( 2*i ))
        end do
     end if

  end do

end subroutine xlarnv
