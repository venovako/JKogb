real(kind=16) function xdot(n, xx, incx, xy, incy)

  implicit none

!     .. Scalar Arguments ..
  integer, intent(in) :: incx, incy, n

!     .. Array Arguments ..
  real(kind=16), intent(in) :: xx(*), xy(*)

!  Purpose
!  =======
!
!     XDOT forms the dot product of two vectors.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack.
!
!  =====================================================================

  real(kind=16), parameter :: zero = +0.0e+0_16

!     .. Local Scalars ..
  real(kind=16) :: xtemp
  integer :: i, ix, iy

  xdot = zero
  xtemp = zero

  if (n .le. 0) return

  ix = 1
  iy = 1
  if (incx .lt. 0) ix = (-n+1)*incx + 1
  if (incy .lt. 0) iy = (-n+1)*incy + 1
  do i = 1, n
     xtemp = xtemp + xx(ix)*xy(iy)
     ix = ix + incx
     iy = iy + incy
  end do
  xdot = xtemp

end function xdot
