subroutine xscal(n, xa, xx, incx)

  implicit none

!     .. Scalar Arguments ..
  real(kind=16), intent(in) :: xa
  integer, intent(in) :: incx, n

!     .. Array Arguments ..
  real(kind=16), intent(inout) :: xx(*)

!  Purpose
!  =======
!
!     XSCAL scales a vector by a constant.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack.
!
!  =====================================================================

!     .. Local Scalars ..
  integer :: i, nincx

  if (n .le. 0 .or. incx .le. 0) return

  nincx = n*incx
  do i = 1, nincx, incx
     xx(i) = xa*xx(i)
  end do

end subroutine xscal
