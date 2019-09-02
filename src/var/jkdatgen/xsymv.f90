subroutine xsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)

  implicit none

!     .. Scalar Arguments ..
  real(kind=16), intent(in) :: alpha, beta
  integer, intent(in) :: incx, incy, lda, n
  character, intent(in) :: uplo

!     .. Array Arguments ..
  real(kind=16), intent(in) :: a(lda,*), x(*)
  real(kind=16), intent(inout) :: y(*)

!  Purpose
!  =======
!
!  XSYMV  performs the matrix-vector  operation
!
!     y := alpha*A*x + beta*y,
!
!  where alpha and beta are scalars, x and y are n element vectors and
!  A is an n by n symmetric matrix.
!
!  Arguments
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the array A is to be referenced as
!           follows:
!
!              UPLO = 'U' or 'u'   Only the upper triangular part of A
!                                  is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the lower triangular part of A
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL(KIND=16).
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - REAL(KIND=16) array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the symmetric matrix and the strictly
!           lower triangular part of A is not referenced.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the symmetric matrix and the strictly
!           upper triangular part of A is not referenced.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!  X      - REAL(KIND=16) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
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
!  Y      - REAL(KIND=16) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y. On exit, Y is overwritten by the updated
!           vector y.
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
  real(kind=16) :: temp1, temp2
  integer :: i, info, ix, iy, j, jx, jy, kx, ky

!     .. External Functions ..
  logical, external :: lsame

!     .. External Subroutines ..
  external :: xerbla

!     Test the input parameters.

  info = 0
  if (.not. lsame(uplo,'U') .and. .not. lsame(uplo,'L')) then
     info = 1
  else if (n .lt. 0) then
     info = 2
  else if (lda .lt. max(1,n)) then
     info = 5
  else if (incx .eq. 0) then
     info = 7
  else if (incy .eq. 0) then
     info = 10
  end if
  if (info .ne. 0) then
     call xerbla('XSYMV ',info)
     return
  end if

!     Quick return if possible.

  if ((n .eq. 0) .or. ((alpha .eq. zero) .and. (beta .eq. one))) return

!     Set up the start points in  X  and  Y.

  if (incx .gt. 0) then
     kx = 1
  else
     kx = 1 - (n-1)*incx
  end if
  if (incy .gt. 0) then
     ky = 1
  else
     ky = 1 - (n-1)*incy
  end if

!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.

!     First form  y := beta*y.

  if (beta .ne. one) then
     if (incy .eq. 1) then
        if (beta .eq. zero) then
           do i = 1, n
              y(i) = zero
           end do
        else
           do i = 1, n
              y(i) = beta*y(i)
           end do
        end if
     else
        iy = ky
        if (beta .eq. zero) then
           do i = 1, n
              y(iy) = zero
              iy = iy + incy
           end do
        else
           do i = 1, n
              y(iy) = beta*y(iy)
              iy = iy + incy
           end do
        end if
     end if
  end if
  if (alpha .eq. zero) return
  if (lsame(uplo,'U')) then

!        Form  y  when A is stored in upper triangle.

     if ((incx .eq. 1) .and. (incy .eq. 1)) then
        do j = 1, n
           temp1 = alpha*x(j)
           temp2 = zero
           do i = 1, j - 1
              y(i) = y(i) + temp1*a(i,j)
              temp2 = temp2 + a(i,j)*x(i)
           end do
           y(j) = y(j) + temp1*a(j,j) + alpha*temp2
        end do
     else
        jx = kx
        jy = ky
        do j = 1, n
           temp1 = alpha*x(jx)
           temp2 = zero
           ix = kx
           iy = ky
           do i = 1, j - 1
              y(iy) = y(iy) + temp1*a(i,j)
              temp2 = temp2 + a(i,j)*x(ix)
              ix = ix + incx
              iy = iy + incy
           end do
           y(jy) = y(jy) + temp1*a(j,j) + alpha*temp2
           jx = jx + incx
           jy = jy + incy
        end do
     end if
  else

!        Form  y  when A is stored in lower triangle.

     if ((incx .eq. 1) .and. (incy .eq. 1)) then
        do j = 1, n
           temp1 = alpha*x(j)
           temp2 = zero
           y(j) = y(j) + temp1*a(j,j)
           do i = j + 1, n
              y(i) = y(i) + temp1*a(i,j)
              temp2 = temp2 + a(i,j)*x(i)
           end do
           y(j) = y(j) + alpha*temp2
        end do
     else
        jx = kx
        jy = ky
        do j = 1, n
           temp1 = alpha*x(jx)
           temp2 = zero
           y(jy) = y(jy) + temp1*a(j,j)
           ix = jx
           iy = jy
           do i = j + 1, n
              ix = ix + incx
              iy = iy + incy
              y(iy) = y(iy) + temp1*a(i,j)
              temp2 = temp2 + a(i,j)*x(ix)
           end do
           y(jy) = y(jy) + alpha*temp2
           jx = jx + incx
           jy = jy + incy
        end do
     end if
  end if

end subroutine xsymv
