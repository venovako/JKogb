subroutine xrot(n, xx, incx, xy, incy, c, s)

  implicit none

!     .. Scalar Arguments ..
  real(kind=16), intent(in) :: c, s
  integer, intent(in) :: incx, incy, n

!     .. Array Arguments ..
  real(kind=16), intent(inout) :: xx(*), xy(*)

!  Purpose
!  =======
!
!     XROT applies a plane rotation.
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
  if (incx .eq. 1 .and. incy .eq. 1) go to 10

!       code for unequal increments or equal increments not equal
!         to 1

  ix = 1
  iy = 1
  if (incx .lt. 0) ix = (-n+1)*incx + 1
  if (incy .lt. 0) iy = (-n+1)*incy + 1
  do i = 1, n
     xtemp = c*xx(ix) + s*xy(iy)
     xy(iy) = c*xy(iy) - s*xx(ix)
     xx(ix) = xtemp
     ix = ix + incx
     iy = iy + incy
  end do
  return

!       code for both increments equal to 1

10 do i = 1, n
     xtemp = c*xx(i) + s*xy(i)
     xy(i) = c*xy(i) - s*xx(i)
     xx(i) = xtemp
  end do

end subroutine xrot
