real(kind=16) function xnrm2(n, x, incx)

  implicit none

!     .. Scalar Arguments ..
  integer, intent(in) :: incx, n

!     .. Array Arguments ..
  real(kind=16), intent(in) :: x(*)

!  Purpose
!  =======
!
!  XNRM2 returns the euclidean norm of a vector via the function
!  name, so that
!
!     XNRM2 := sqrt( x'*x )
!
!  Further Details
!  ===============
!
!     Sven Hammarling, Nag Ltd.
!
!  =====================================================================

!     .. Parameters ..
  real(kind=16), parameter :: one = +1.0e+0_16, zero = +0.0e+0_16

!     .. Local Scalars ..
  real(kind=16) :: absxi, norm, scal, ssq
  integer :: ix

  if (n .lt. 1 .or. incx .lt. 1) then
     norm = zero
  else if (n .eq. 1) then
     norm = abs(x(1))
  else
     scal = zero
     ssq = one

!        The following loop is equivalent to this call to the LAPACK
!        auxiliary routine:
!        CALL XLASSQ( N, X, INCX, SCALE, SSQ )

     do ix = 1, 1 + (n-1)*incx, incx
        if (x(ix) .ne. zero) then
           absxi = abs(x(ix))
           if (scal .lt. absxi) then
              ssq = one + ssq * (scal/absxi)**2
              scal = absxi
           else
              ssq = ssq + (absxi/scal)**2
           end if
        end if
     end do
     norm = scal * sqrt(ssq)
  end if

  xnrm2 = norm

end function xnrm2
