subroutine xcopy(n, xx, incx, xy, incy)

  implicit none

!     .. Scalar Arguments ..
  integer, intent(in) :: incx, incy, n

!     .. Array Arguments ..
  real(kind=16), intent(in) :: xx(*)
  real(kind=16), intent(out) :: xy(*)

!  Purpose
!  =======
!
!     XCOPY copies a vector, x, to a vector, y.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack.
!
!  =====================================================================

!     .. Local Scalars ..
  integer :: i, ix, iy

  if (n .le. 0) return

!        code for unequal increments or equal increments
!          not equal to 1

  ix = 1
  iy = 1
  if (incx .lt. 0) ix = (-n+1)*incx + 1
  if (incy .lt. 0) iy = (-n+1)*incy + 1
  do i = 1, n
     xy(iy) = xx(ix)
     ix = ix + incx
     iy = iy + incy
  end do

end subroutine xcopy
