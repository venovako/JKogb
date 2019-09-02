integer function ixamax(n, xx, incx)

  implicit none

!     .. Scalar Arguments ..
  integer, intent(in) :: incx, n

!     .. Array Arguments ..
  real(kind=16), intent(in) :: xx(*)

!  Purpose
!  =======
!
!     IXAMAX finds the index of element having max. absolute value.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack.
!
!  =====================================================================

!     .. Local Scalars ..
  real(kind=16) :: xmax
  integer :: i, ix

  ixamax = 0
  if (n .lt. 1 .or. incx .le. 0) return
  ixamax = 1
  if (n .eq. 1) return
  if (incx .eq. 1) go to 10

!        code for increment not equal to 1

  ix = 1
  xmax = abs(xx(1))
  ix = ix + incx
  do i = 2, n
     if (abs(xx(ix)) .le. xmax) go to 5
     ixamax = i
     xmax = abs(xx(ix))
5    ix = ix + incx
  end do
  return

!        code for increment equal to 1

10 xmax = abs(xx(1))
  do i = 2, n
     if (abs(xx(i)) .le. xmax) cycle
     ixamax = i
     xmax = abs(xx(i))
  end do

end function ixamax
