subroutine xaxpy(n, xa, xx, incx, xy, incy)

  implicit none

!     .. Scalar Arguments ..
  real(kind=16), intent(in) :: xa
  integer, intent(in) :: incx, incy, n

!     .. Array Arguments ..
  real(kind=16), intent(in) :: xx(*)
  real(kind=16), intent(inout) :: xy(*)

!  Purpose
!  =======
!
!     XAXPY constant times a vector plus a vector.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack.
!
!  =====================================================================

  real(kind=16), parameter :: zero = +0.0e+0_16

!     .. Local Scalars ..
  integer :: i, ix, iy

  if (n .le. 0) return
  if (xa .eq. zero) return

  ix = 1
  iy = 1
  if (incx .lt. 0) ix = (-n+1)*incx + 1
  if (incy .lt. 0) iy = (-n+1)*incy + 1
  do i = 1, n
     xy(iy) = xy(iy) + xa*xx(ix)
     ix = ix + incx
     iy = iy + incy
  end do

end subroutine xaxpy
