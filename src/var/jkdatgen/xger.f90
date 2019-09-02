subroutine xger(m, n, alpha, x, incx, y, incy, a, lda)

  implicit none

!     .. Scalar Arguments ..
  real(kind=16), intent(in) :: alpha
  integer, intent(in) :: incx, incy, lda, m, n

!     .. Array Arguments ..
  real(kind=16), intent(inout) :: a(lda,*)
  real(kind=16), intent(in) :: x(*), y(*)

!  Purpose
!  =======
!
!  XGER   performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Arguments
!  ==========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL(KIND=16).
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - REAL(KIND=16) array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - REAL(KIND=16) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - REAL(KIND=16) array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  Further Details
!  ===============
!
!  Level 2 Blas routine.
!
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!  =====================================================================

!     .. Parameters ..
  real(kind=16), parameter :: zero = +0.0e+0_16

!     .. Local Scalars ..
  real(kind=16) :: temp
  integer :: i, info, ix, j, jy, kx

!     .. External Subroutines ..
  external :: xerbla

!     Test the input parameters.

  info = 0
  if (m .lt. 0) then
     info = 1
  else if (n .lt. 0) then
     info = 2
  else if (incx .eq. 0) then
     info = 5
  else if (incy .eq. 0) then
     info = 7
  else if (lda .lt. max(1,m)) then
     info = 9
  end if
  if (info .ne. 0) then
     call xerbla('XGER  ',info)
     return
  end if

!     Quick return if possible.

  if ((m .eq. 0) .or. (n .eq. 0) .or. (alpha .eq. zero)) return

!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.

  if (incy .gt. 0) then
     jy = 1
  else
     jy = 1 - (n-1)*incy
  end if
  if (incx .eq. 1) then
     do j = 1, n
        if (y(jy) .ne. zero) then
           temp = alpha*y(jy)
           do i = 1, m
              a(i,j) = a(i,j) + x(i)*temp
           end do
        end if
        jy = jy + incy
     end do
  else
     if (incx .gt. 0) then
        kx = 1
     else
        kx = 1 - (m-1)*incx
     end if
     do j = 1, n
        if (y(jy) .ne. zero) then
           temp = alpha*y(jy)
           ix = kx
           do i = 1, m
              a(i,j) = a(i,j) + x(ix)*temp
              ix = ix + incx
           end do
        end if
        jy = jy + incy
     end do
  end if

end subroutine xger
