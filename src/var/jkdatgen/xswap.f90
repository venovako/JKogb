subroutine xswap(n, xx, incx, xy, incy)

  implicit none

!     .. Scalar Arguments ..
  integer, intent(in) :: incx, incy, n

!     .. Array Arguments ..
  real(kind=16), intent(inout) :: xx(*), xy(*)

!  Purpose
!  =======
!
!     interchanges two vectors.
!     uses unrolled loops for increments equal one.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack.
!
!  =====================================================================

!     .. Local Scalars ..
  real(kind=16) :: xtemp
  integer :: i, ix, iy

  if (n .le. 0) return

  ix = 1
  iy = 1
  if (incx .lt. 0) ix = (-n+1)*incx + 1
  if (incy .lt. 0) iy = (-n+1)*incy + 1
  do i = 1, n
     xtemp = xx(ix)
     xx(ix) = xy(iy)
     xy(iy) = xtemp
     ix = ix + incx
     iy = iy + incy
  end do

end subroutine xswap
