subroutine xsyr(uplo, n, alpha, x, incx, a, lda)

  implicit none

!     .. Scalar Arguments ..
  real(kind=16), intent(in) :: alpha
  integer, intent(in) :: incx, lda, n
  character, intent(in) :: uplo

!     .. Array Arguments ..
  real(kind=16), intent(inout) :: a(lda,*)
  real(kind=16), intent(in) :: x(*)

!  Purpose
!  =======
!
!  XSYR   performs the symmetric rank 1 operation
!
!     A := alpha*x*x' + A,
!
!  where alpha is a real scalar, x is an n element vector and A is an
!  n by n symmetric matrix.
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
  real(kind=16) :: temp
  integer :: i, info, ix, j, jx, kx

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
  else if (lda .lt. max(1,n)) then
     info = 7
  end if
  if (info .ne. 0) then
     call xerbla('XSYR  ',info)
     return
  end if

!     Quick return if possible.

  if ((n .eq. 0) .or. (alpha .eq. zero)) return

!     Set the start point in X if the increment is not unity.

  if (incx .le. 0) then
     kx = 1 - (n-1)*incx
  else if (incx .ne. 1) then
     kx = 1
  end if

!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.

  if (lsame(uplo,'U')) then

!        Form  A  when A is stored in upper triangle.

     if (incx .eq. 1) then
        do j = 1, n
           if (x(j) .ne. zero) then
              temp = alpha * x(j)
              do i = 1, j
                 a(i,j) = a(i,j) + x(i) * temp
              end do
           end if
        end do
     else
        jx = kx
        do j = 1, n
           if (x(jx) .ne. zero) then
              temp = alpha * x(jx)
              ix = kx
              do i = 1, j
                 a(i,j) = a(i,j) + x(ix) * temp
                 ix = ix + incx
              end do
           end if
           jx = jx + incx
        end do
     end if
  else

!        Form  A  when A is stored in lower triangle.

     if (incx .eq. 1) then
        do j = 1, n
           if (x(j) .ne. zero) then
              temp = alpha * x(j)
              do i = j, n
                 a(i,j) = a(i,j) + x(i) * temp
              end do
           end if
        end do
     else
        jx = kx
        do j = 1, n
           if (x(jx) .ne. zero) then
              temp = alpha * x(jx)
              ix = jx
              do i = j, n
                 a(i,j) = a(i,j) + x(ix) * temp
                 ix = ix + incx
              end do
           end if
           jx = jx + incx
        end do
     end if
  end if

end subroutine xsyr
