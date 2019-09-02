subroutine xgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)

  implicit none

!     .. Scalar Arguments ..
  real(kind=16), intent(in) :: alpha, beta
  integer, intent(in) :: incx, incy, lda, m, n
  character, intent(in) :: trans

!     .. Array Arguments ..
  real(kind=16), intent(in) :: a(lda,*), x(*)
  real(kind=16), intent(inout) :: y(*)

!  Purpose
!  =======
!
!  XGEMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
!  Arguments
!  ==========
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
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
!  A      - REAL(KIND=16) array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - REAL(KIND=16) array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - REAL(KIND=16).
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - REAL(KIND=16) array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
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
  real(kind=16), parameter :: one = +1.0e+0_16, zero = +0.0e+0_16

!     .. Local Scalars ..
  real(kind=16) :: temp
  integer :: i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny

!     .. External Functions ..
  logical, external :: lsame

!     .. External Subroutines ..
  external :: xerbla

!     Test the input parameters.

  info = 0

  if (.not. lsame(trans,'N') .and. .not. lsame(trans,'T') .and. .not. lsame(trans,'C')) then
     info = 1
  else if (m .lt. 0) then
     info = 2
  else if (n .lt. 0) then
     info = 3
  else if (lda .lt. max(1,m)) then
     info = 6
  else if (incx .eq. 0) then
     info = 8
  else if (incy .eq. 0) then
     info = 11
  end if
  if (info .ne. 0) then
     call xerbla('XGEMV ',info)
     return
  end if

!     Quick return if possible.

  if ((m .eq. 0) .or. (n .eq. 0) .or. ((alpha .eq. zero) .and. (beta .eq. one))) return

!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.

  if (lsame(trans,'N')) then
     lenx = n
     leny = m
  else
     lenx = m
     leny = n
  end if
  if (incx .gt. 0) then
     kx = 1
  else
     kx = 1 - (lenx-1)*incx
  end if
  if (incy .gt. 0) then
     ky = 1
  else
     ky = 1 - (leny-1)*incy
  end if

!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.

!     First form  y := beta*y.

  if (beta .ne. one) then
     if (incy .eq. 1) then
        if (beta .eq. zero) then
           do i = 1, leny
              y(i) = zero
           end do
        else
           do i = 1, leny
              y(i) = beta*y(i)
           end do
        end if
     else
        iy = ky
        if (beta .eq. zero) then
           do i = 1, leny
              y(iy) = zero
              iy = iy + incy
           end do
        else
           do i = 1, leny
              y(iy) = beta*y(iy)
              iy = iy + incy
           end do
        end if
     end if
  end if
  if (alpha .eq. zero) return
  if (lsame(trans,'N')) then

!        Form  y := alpha*A*x + y.

     jx = kx
     if (incy .eq. 1) then
        do j = 1, n
           if (x(jx) .ne. zero) then
              temp = alpha*x(jx)
              do i = 1, m
                 y(i) = y(i) + temp*a(i,j)
              end do
           end if
           jx = jx + incx
        end do
     else
        do j = 1, n
           if (x(jx) .ne. zero) then
              temp = alpha*x(jx)
              iy = ky
              do i = 1, m
                 y(iy) = y(iy) + temp*a(i,j)
                 iy = iy + incy
              end do
           end if
           jx = jx + incx
        end do
     end if
  else

!        Form  y := alpha*A'*x + y.

     jy = ky
     if (incx .eq. 1) then
        do j = 1, n
           temp = zero
           do i = 1, m
              temp = temp + a(i,j)*x(i)
           end do
           y(jy) = y(jy) + alpha*temp
           jy = jy + incy
        end do
     else
        do j = 1, n
           temp = zero
           ix = kx
           do i = 1, m
              temp = temp + a(i,j)*x(ix)
              ix = ix + incx
           end do
           y(jy) = y(jy) + alpha*temp
           jy = jy + incy
        end do
     end if
  end if

end subroutine xgemv
