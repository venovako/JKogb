subroutine xsybpc( uplo, n, a, lda, nrank, diagj, ipiv, jvec, info )

  !     Modified by Singers, January 8, 2006.
  !
  !     Purpose
  !     =======
  !
  !     XSYBPC computes the modified Bunch-Parlett factorization of a N-by-N
  !     real symmetric matrix A
  !
  !     A = G * J * transpose( G ).
  !
  !     Arguments
  !     =========
  !
  !     UPLO    (input) CHARACTER*1
  !     = 'U':  Upper triangle of A is stored;
  !     = 'L':  Lower triangle of A is stored.
  !
  !     N       (input) INTEGER
  !     Order of the input matrix A, N >= 0.
  !
  !     A       (input/output) REAL(KIND=16) array, dimension (LDA,N)
  !     On entry, the symmetric matrix A.  If UPLO = 'U', the leading
  !     N-by-N upper triangular part of A contains the upper
  !     triangular part of the matrix A. If UPLO = 'L', the
  !     leading N-by-N lower triangular part of A contains the lower
  !     triangular part of the matrix A.
  !
  !     On exit, A contains the matrix G.
  !
  !     LDA     (input) INTEGER
  !     The leading dimension of the array A.  LDA >= max(1,N).
  !
  !     NRANK   (output) INTEGER
  !     Contains the rank of A.
  !
  !     DIAGJ   (workspace) REAL(KIND=16) array, dimension(N)
  !     Contains the signs of the eigenvalues.
  !
  !     IPIV    (workspace) INTEGER array, dimension(N)
  !
  !     JVEC    (output) INTEGER array, dimension(N)
  !     Contains the diagonal of the matrix J; JVEC( I ) = 1 or -1.
  !     If NRANK < N, only the first NRANK values are set.
  !
  !     INFO    (output) INTEGER
  !     = 0:  successful exit - the first N columns of the array A
  !     contain the matrix G (full column rank).
  !     < 0:  if INFO = -i, with 0 < i <= 1000, the i-th argument
  !     had an illegal value;
  !     > 0:  some eigenvalues are zero and INFO specifies the
  !     rank of A.

  implicit none

  character(len=1), intent(in) :: uplo
  integer, intent(in) :: n
  real(kind=16), dimension(lda, n), intent(inout) :: a
  integer, intent(in) :: lda
  integer, intent(out) :: nrank
  real(kind=16), dimension(n), intent(out) :: diagj
  integer, dimension(n), intent(out) :: ipiv
  integer, dimension(n), intent(out) :: jvec
  integer, intent(out) :: info

  real(kind=16), parameter :: zero = +0.0e+0_16

  logical :: upper
  integer :: i, infod
  integer :: k, kp, nr2

  logical, external :: lsame
  external :: xcopy, xlaset, xswap, xsyjf2

  !     Test the input arguments
  info = 0
  upper = lsame( uplo, 'U' )

  if ( .not. upper .and. .not. lsame( uplo, 'L' ) ) then
     info = -1
  else if ( n .lt. 0 ) then
     info = -2
  else if ( lda .lt. max( 1, n ) ) then
     info = -4
  end if

  if ( info .ne. 0 ) return

  !     Quick return, if possible.
  if ( n .eq. 0 ) then
     nrank = 0
     return
  end if

  !     Compute the factorization from the lower triangle of A.
  if ( upper ) then
     !     Copy the upper triangle to the lower triangle
     do i = 1, n - 1
        call xcopy( n-i, a( i, i+1 ), lda, a( i+1, i ), 1 )
     end do
  end if

  !     Compute the factorization A = L*J*L', where L is a product
  !     of permutation and lower block triangular matrices with
  !     1-by-1 and 2-by-2 diagonal blocks, L' is the transpose of L,
  !     and J is diagonal with diagonal elements equal to 1 or -1.
  call xsyjf2( 'L', n, a, lda, diagj, ipiv, infod )

  !     Set NRANK to the rank of A
  if ( infod .eq. 0 ) then
     nrank = n
  else
     nrank = infod - 1
  end if

  !     Create the factor G such that  P*A*P' = G*J*G', where G is
  !     lower block triangular matrix, and J is as above.

  !     Empty the the the array A above the second superdiagonal
  nr2 = nrank - 2
  if ( nr2 .gt. 0 ) call xlaset( 'U', nr2, nr2, zero, zero, a( 1, 3 ), lda )

  !     Set the element A( 1, 2 ) or A( 2, 3 ) to zero if necessary
  if ( ipiv( 1 ) .gt. 0 ) then
     if ( nrank .gt. 1 ) a( 1, 2 ) = zero
     !     K is the loop index, increasing from 2 or 3 to NRANK in steps
     !     of 1 or 2, depending on the size of the diagonal blocks.
     k = 2
  else
     if ( nrank .gt. 2 ) a( 2, 3 ) = zero
     k = 3
  end if

20 continue
  !     If K > NRANK, exit from loop.
  if ( k .gt. nrank ) go to 30

  if ( ipiv( k ) .gt. 0 ) then
     !     1 x 1 diagonal block
     !
     !     Interchange rows K and IPIV(K).
     kp = ipiv( k )
     if ( kp .ne. k ) call xswap( k-1, a( k, 1 ), lda, a( kp, 1 ), lda )

     !     Set the element A( k, k+1 ) of the upper triangle to zero
     if ( nrank .gt. k ) a( k, k+1 ) = zero
     k = k + 1
  else
     !     2 x 2 diagonal block
     !
     !     Interchange rows K and -IPIV(K).
     kp = -ipiv( k )
     if ( kp .ne. k ) call xswap( k-1, a( k, 1 ), lda, a( kp, 1 ), lda )

     !     Interchange rows K+1 and -IPIV(K+1).
     kp = -ipiv( k+1 )
     if ( kp .ne. k+1 ) call xswap( k-1, a( k+1, 1 ), lda, a( kp, 1 ), lda )

     !     Set the element A( K+1, K+2 ) of the upper triangle to zero
     if ( nrank .gt. k+1 ) a( k+1, k+2 ) = zero
     k = k + 2
  end if
  go to 20

30 continue
  do i = 1, nrank
     jvec( i ) = int( diagj( i ) )
  end do

  !     If column rank defect occured, set INFO = RANK
  if ( nrank .lt. n ) info = nrank

end subroutine xsybpc
