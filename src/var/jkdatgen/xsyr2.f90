subroutine xsyr2(uplo, n, alpha, x, incx, y, incy, a, lda)

  implicit none

!     .. Scalar Arguments ..
  real(kind=16), intent(in) :: alpha
  integer, intent(in) :: incx, incy, lda, n
  character, intent(in) :: uplo

!    .. Array Arguments ..
  real(kind=16), intent(inout) :: a(lda,*)
  real(kind=16), intent(in) :: x(*), y(*)

!  Purpose
!  =======
!
!  XSYR2  performs the symmetric rank 2 operation
!
!     A := alpha*x*y' + alpha*y*x' + A,
!
!  where alpha is a scalar, x and y are n element vectors and A is an n
!  by n symmetric matrix.
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
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the symmetric matrix and the strictly
!           lower triangular part of A is not referenced. On exit, the
!           upper triangular part of the array A is overwritten by the
!           upper triangular part of the updated matrix.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the symmetric matrix and the strictly
!           upper triangular part of A is not referenced. On exit, the
!           lower triangular part of the array A is overwritten by the
!           lower triangular part of the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
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
  else if (incx .eq. 0) then
     info = 5
  else if (incy .eq. 0) then
     info = 7
  else if (lda .lt. max(1,n)) then
     info = 9
  end if
  if (info .ne. 0) then
     call xerbla('XSYR2 ',info)
     return
  end if

!     Quick return if possible.

  if ((n .eq. 0) .or. (alpha .eq. zero)) return

!     Set up the start points in X and Y if the increments are not both
!     unity.

  if ((incx .ne. 1) .or. (incy .ne. 1)) then
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
     jx = kx
     jy = ky
  end if

!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.

  if (lsame(uplo,'U')) then

!        Form  A  when A is stored in the upper triangle.

     if ((incx .eq. 1) .and. (incy .eq. 1)) then
        do j = 1, n
           if ((x(j).ne.zero) .or. (y(j).ne.zero)) then
              temp1 = alpha*y(j)
              temp2 = alpha*x(j)
              do i = 1, j
                 a(i,j) = a(i,j) + x(i)*temp1 + y(i)*temp2
              end do
           end if
        end do
     else
        do j = 1, n
           if ((x(jx) .ne. zero) .or. (y(jy) .ne. zero)) then
              temp1 = alpha*y(jy)
              temp2 = alpha*x(jx)
              ix = kx
              iy = ky
              do i = 1, j
                 a(i,j) = a(i,j) + x(ix)*temp1 + y(iy)*temp2
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           jx = jx + incx
           jy = jy + incy
        end do
     end if
  else
     
!        Form  A  when A is stored in the lower triangle.
     
     if ((incx .eq. 1) .and. (incy .eq. 1)) then
        do j = 1, n
           if ((x(j) .ne. zero) .or. (y(j) .ne. zero)) then
              temp1 = alpha*y(j)
              temp2 = alpha*x(j)
              do i = j, n
                 a(i,j) = a(i,j) + x(i)*temp1 + y(i)*temp2
              end do
           end if
        end do
     else
        do j = 1, n
           if ((x(jx) .ne. zero) .or. (y(jy) .ne. zero)) then
              temp1 = alpha*y(jy)
              temp2 = alpha*x(jx)
              ix = jx
              iy = jy
              do i = j, n
                 a(i,j) = a(i,j) + x(ix)*temp1 + y(iy)*temp2
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           jx = jx + incx
           jy = jy + incy
        end do
     end if
  end if

end subroutine xsyr2
