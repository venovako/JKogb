real(kind=8) function dlaran( iseed )

  implicit none

!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..

!     .. Array Arguments ..
  integer, intent(inout) :: iseed(4)

!  Purpose
!  =======
!
!  DLARAN returns a random real number from a uniform (0,1)
!  distribution.
!
!  Arguments
!  =========
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
!  This routine uses a multiplicative congruential method with modulus
!  2**48 and multiplier 33952834046453 (see G.S.Fishman,
!  'Multiplicative congruential random number generators with modulus
!  2**b: an exhaustive analysis for b = 32 and a partial analysis for
!  b = 48', Math. Comp. 189, pp 331-344, 1990).
!
!  48-bit integers are stored in 4 integer array elements with 12 bits
!  per element. Hence the routine is portable across machines with
!  integers of 32 bits or more.
!
!  =====================================================================

!     .. Parameters ..
  integer, parameter :: m1 = 494, m2 = 322, m3 = 2508, m4 = 2549
  real(kind=8), parameter :: one = +1.0e+0_8
  integer, parameter :: ipw2 = 4096
  real(kind=8), parameter :: r = one / ipw2

!     .. Local Scalars ..
  integer :: it1, it2, it3, it4
  real(kind=8) :: rndout

!     .. Executable Statements ..
10 continue

!     multiply the seed by the multiplier modulo 2**48

  it4 = iseed( 4 )*m4
  it3 = it4 / ipw2
  it4 = it4 - ipw2*it3
  it3 = it3 + iseed( 3 )*m4 + iseed( 4 )*m3
  it2 = it3 / ipw2
  it3 = it3 - ipw2*it2
  it2 = it2 + iseed( 2 )*m4 + iseed( 3 )*m3 + iseed( 4 )*m2
  it1 = it2 / ipw2
  it2 = it2 - ipw2*it1
  it1 = it1 + iseed( 1 )*m4 + iseed( 2 )*m3 + iseed( 3 )*m2 + iseed( 4 )*m1
  it1 = mod( it1, ipw2 )

!     return updated seed

  iseed( 1 ) = it1
  iseed( 2 ) = it2
  iseed( 3 ) = it3
  iseed( 4 ) = it4

!     convert 48-bit integer to a real number in the interval (0,1)

  rndout = r*(real(it1,8) + r*(real(it2,8) + r*(real(it3,8) + r*(real(it4,8)))))

  if (rndout .eq. one) then
!        If a real number has n bits of precision, and the first
!        n bits of the 48-bit integer above happen to be all 1 (which
!        will occur about once every 2**n calls), then DLARAN will
!        be rounded to exactly 1.0. 
!        Since DLARAN is not supposed to return exactly 0.0 or 1.0
!        (and some callers of DLARAN, such as CLARND, depend on that),
!        the statistically correct thing to do in this situation is
!        simply to iterate again.
!        N.B. the case DLARAN = 0.0 should not be possible.

     goto 10
  end if

  dlaran = rndout

end function dlaran
