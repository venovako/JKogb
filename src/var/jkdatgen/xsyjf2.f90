subroutine xsyjf2( uplo, n, a, lda, diagj, ipiv, info )

  !     Modified by Singers, January 8, 2006.
  !
  !     Proposal by Ivan Slapnicar
  !     University of Split, Croatia
  !     slap@split.fesb.hr
  !     December 15, 1993
  !
  !     Purpose
  !     =======
  !
  !     XSYJF2 computes the factorization of a real symmetric matrix A using
  !     a modification of the Bunch-Kaufman diagonal pivoting method:
  !
  !     A = U*J*U'  or  A = L*J*L'
  !
  !     where U (or L) is a product of permutation and  upper (lower) block
  !     triangular matrices with 1-by-1 and 2-by-2 diagonal blocks, U' is
  !     the transpose of U, and J is diagonal with diagonal elements equal
  !     to 1 or -1 (see below for further details).
  !
  !     This is the unblocked version of the algorithm, calling Level 2 BLAS.
  !
  !     Arguments
  !     =========
  !
  !     UPLO    (input) CHARACTER*1
  !     Specifies whether the upper or lower triangular part of the
  !     symmetric matrix A is stored:
  !     = 'U':  Upper triangular
  !     = 'L':  Lower triangular
  !
  !     N       (input) INTEGER
  !     The order of the matrix A.  N >= 0.
  !
  !     A       (input/output) REAL(KIND=16) array, dimension (LDA,N)
  !     On entry, the symmetric matrix A.  If UPLO = 'U', the leading
  !     N-by-N upper triangular part of A contains the upper
  !     triangular part of the matrix A. If UPLO = 'L', the
  !     leading N-by-N lower triangular part of A contains the lower
  !     triangular part of the matrix A.
  !
  !     On exit, the multipliers used to obtain the block triangular
  !     factor U or L (see below for further details).
  !
  !     LDA     (input) INTEGER
  !     The leading dimension of the array A.  LDA >= max(1,N).
  !
  !     DIAGJ   (output) REAL(KIND=16) array, dimension(N)
  !     The diagonal of the matrix J.
  !
  !     IPIV    (output) INTEGER array, dimension (N)
  !     Details of the interchanges and the block structure of D.
  !     If IPIV(k) > 0, then rows and columns k and IPIV(k) were
  !     interchanged and in the k-th step a 1-by-1 diagonal block
  !     was used (see below for further details).
  !     If UPLO = 'U' and IPIV(k) < 0 and IPIV(k-1) < 0, then rows
  !     and columns k and -IPIV(k) and k-1 and -IPIV(k-1) were
  !     interchanged, and in the k-th step a 2-by-2 diagonal block
  !     was used.
  !     If UPLO = 'L' and IPIV(k) < 0 and IPIV(k+1) < 0, then rows
  !     and columns k and -IPIV(k) and k+1 and -IPIV(k+1) were
  !     interchanged, and in the k-th step a 2-by-2 diagonal block
  !     was used.
  !
  !     INFO    (output) INTEGER
  !     = 0: successful exit
  !     < 0: if INFO = -k, the k-th argument had an illegal value
  !     > 0: if UPLO = 'U' and INFO = k, then in the k-th step the
  !     leading submatrix A(1:k,1:k) was exactly zero, and the
  !     rank of A equals n-k;
  !     if UPLO = 'L' and INFO = k, then in the k-th step the
  !     trailing submatrix A(k:n,k:n) was exactly zero, and the
  !     rank of A equals k-1;
  !
  !     Further Details
  !     ===============
  !
  !     If UPLO = 'U', then A = U*J*U', where
  !     U = P(n)*U(n)* ... *P(k)U(k)* ...,
  !     i.e., U is a product of terms P(k)*U(k), where k decreases from n to
  !     1 in steps of 1 or 2, and J is a diagonal matrix with diagonal
  !     elements equal to 1 or -1. P(k) is a permutation matrix as defined
  !     by IPIV(k), and U(k) is a upper block triangular matrix, such that
  !
  !             (   I    v    0   )   k-s
  !     U(k) =  (   0    d    0   )   s
  !             (   0    0    I   )   n-k
  !                k-s   s   n-k
  !
  !     Here s = 1 or 2, and ( v ) overwrites A(1:k,k-s+1:k).
  !                          ( d )
  !     
  !     If UPLO = 'L', then A = L*J*L', where
  !     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
  !     i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
  !     n in steps of 1 or 2, and J is a diagonal matrix with diagonal
  !     elements equal to 1 or -1. P(k) is a permutation matrix as defined
  !     by IPIV(k), and L(k) is a lower block triangular matrix, such that
  !
  !             (   I    0     0   )  k-1
  !     L(k) =  (   0    d     0   )  s
  !             (   0    v     I   )  n-k-s+1
  !                k-1   s  n-k-s+1
  !
  !     Here s = 1 or 2, and ( d ) overwrites A(k:n,k:k+s-1).
  !                          ( v )
  !
  !
  !     The following comment is related to a comment in the eigenreduction
  !     routine SSYEJA: the output format
  !
  !     U = P(n)*U(n)* ... *P(k)U(k)* ...,  or
  !     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
  !
  !     is taken from the Lapack routine SSYTF2. For our purpose it might be
  !     more convenient to have the output in the form
  !     P*H*P' = U*J*U'   or  P*H*P' = L*J*L'
  !     as in the original Bunch and Parlett paper from Wilkinson, Reinsch book,
  !     and in the former routine SGJGT. Then the first backward permutation
  !     in SSYEJA would not be necessary. Thus, the final choice is an open
  !     question.

  implicit none

  character(len=1), intent(in) :: uplo
  integer, intent(in) :: n
  real(kind=16), dimension(lda, n), intent(inout) :: a
  integer, intent(in) :: lda
  real(kind=16), dimension(n), intent(out) :: diagj
  integer, dimension(n), intent(out) :: ipiv
  integer, intent(out) :: info

  real(kind=16), parameter ::   zero =  +0.0e+0_16
  real(kind=16), parameter ::    one =  +1.0e+0_16
  real(kind=16), parameter ::  eight =  +8.0e+0_16
  real(kind=16), parameter :: sevten = +17.0e+0_16

  logical :: upper
  integer :: i, j, imax, jmax, k, kdiag, kp, kstep
  real(kind=16) :: alpha, c, diamax, offmax, r1, r2, s, t, temp

  integer, external :: ixamax
  logical, external :: lsame

  external :: xlaev2, xrot, xscal, xswap, xsyr

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

  !     Initialize ALPHA for use in choosing pivot block size.
  alpha = ( one + sqrt( sevten ) ) / eight
  imax = 0
  jmax = 0

  if ( upper ) then
     !     Factorize A as U*J*U' using the upper triangle of A
     !
     !     K is the main loop index, decreasing from N to 1 in steps of 1 or 2
     k = n
10   continue
     !     If K < 1, exit from loop
     if ( k .lt. 1 ) go to 70

     !     Determine rows and columns to be interchanged and whether
     !     a 1-by-1 or 2-by-2 pivot block will be used

     !     KDIAG is the index of the largest diagonal element, and
     !     DIAMAX is its absolute value
     kdiag = ixamax( k, a( 1, 1 ), lda + 1 )
     diamax = abs( a( kdiag, kdiag ) )

     !     IMAX and JMAX are the row- and column-indices of the largest
     !     off-diagonal element, and OFFMAX is its absolute value
     offmax = zero
     if ( k .gt. 1 ) then
        do j = 2, k
           do i = 1, j - 1
              temp = abs( a( i, j ) )
              if ( temp .gt. offmax ) then
                 offmax = temp
                 imax = i
                 jmax = j
              end if
           end do
        end do
     end if

     if ( max( diamax, offmax ) .eq. zero ) then
        !     The rest of the matrix is zero: set INFO and return
        info = k
        go to 70
     end if

     if ( diamax .ge. alpha * offmax ) then
        !     Use 1-by-1 pivot block
        kstep = 1
        kp = kdiag
     else
        !     Use 2-by-2 pivot block
        kstep = 2
        kp = jmax
     end if

     if ( kp .ne. k ) then
        !     Interchange rows and columns K and KP in the leading
        !     submatrix A(1:k,1:k)
        call xswap( kp-1, a( 1, k ), 1, a( 1, kp ), 1 )
        call xswap( k-kp-1, a( kp+1, k ), 1, a( kp, kp+1 ), lda )
        t = a( k, k )
        a( k, k ) = a( kp, kp )
        a( kp, kp ) = t
     end if

     if ( ( kstep .eq. 2 ) .and. ( imax .ne. k-1 ) ) then
        !     Interchange rows and columns K-1 and IMAX in the leading
        !     submatrix A(1:k,1:k)
        call xswap( imax-1, a( 1, k-1 ), 1, a( 1, imax ), 1 )
        call xswap( k-imax-2, a( imax+1, k-1 ), 1, a( imax, imax+1 ), lda )
        t = a( k-1, k-1 )
        a( k-1, k-1 ) = a( imax, imax )
        a( imax, imax ) = t
        t = a( imax, k )
        a( imax, k ) = a( k-1, k )
        a( k-1, k ) = t
     end if

     !     Update the leading submatrix
     if ( kstep .eq. 1 ) then
        !     1-by-1 pivot block D(k): column k now holds
        !
        !     W(1:k-1,k) = U(1:k-1,k) * sqrt( abs(D(k)) ),
        !     W(k,k) = U(k,k) * sqrt( abs(D(k)) ) * J(k),
        !
        !     where U(1:k,k) is the k-th column of U, and
        !     J(k) = sign( D(k) )
        !
        !     Perform a rank-1 update of A(1:k-1,1:k-1) as
        !
        !     A := A - U(k)*J(k,k)*U(k)' = A - W(k)*1/D(k)*W(k)'
        r1 = one / a( k, k )
        if ( k .gt. 1 ) call xsyr( uplo, k-1, -r1, a( 1, k ), 1, a, lda )

        !     Compute the k-th diagonal element of the matrix J
        if ( r1 .ge. zero ) then
           diagj( k ) = one
        else
           diagj( k ) = -one
        end if

        !     Store U(k) in column k
        r1 = sqrt( abs( r1 ) )
        call xscal( k-1, r1, a( 1, k ), 1 )
        a( k, k ) = diagj( k ) / r1
     else
        !     2-by-2 pivot block D(k): let
        !
        !     D(k) = Q(k)**T * X(k) * Q(k)
        !
        !     be the eigendecomposition of D(k), X(k) = diag(R1,R2).
        !     Columns k and k-1 now hold
        !
        !     ( W(1:k-2,k-1) W(1:k-2,k) ) =
        !
        !     ( U(1:k-2,k-1) U(1:k-2,k) )*sqrt(abs(X(k)))*Q(k)**T,
        !
        !     W(k-1:k,k-1:k) =
        !
        !     U(k-1:k,k-1:k)*Q(k)*inv(sqrt(abs(X(k))))*J(k),
        !
        !     where U(k) and U(k-1) are the k-th and (k-1)-th columns
        !     of U, and J(k) = diag( sign(R1), sign(R2) ).
        !
        !     Perform a rank-2 update of A(1:k-2,1:k-2) as
        !
        !     A := A - ( U(k-1) U(k) )*J(k)*( U(k-1) U(k) )'
        !     = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )'
        !
        !     Convert this to two rank-1 updates by using the eigen-
        !     decomposition of D(k)
        call xlaev2( a( k-1, k-1 ), a( k-1, k ), a( k, k ), r1, r2, c, s )
        r1 = one / r1
        r2 = one / r2
        call xrot( k-2, a( 1, k-1 ), 1, a( 1, k ), 1, c, s )
        call xsyr( uplo, k-2, -r1, a( 1, k-1 ), 1, a, lda )
        call xsyr( uplo, k-2, -r2, a( 1, k ), 1, a, lda )

        !     Compute the (k-1)-st and k-th diagonal element of the matrix J
        if ( r1 .ge. zero ) then
           diagj( k-1 ) = one
           diagj( k ) = -one
        else
           diagj( k-1 ) = -one
           diagj( k ) = one
        end if

        !     Store U(k) and U(k-1) in columns k and k-1
        r1 = sqrt( abs( r1 ) )
        r2 = sqrt( abs( r2 ) )
        call xscal( k-2, r1, a( 1, k-1 ), 1 )
        call xscal( k-2, r2, a( 1, k ), 1 )

        r1 = diagj( k-1 ) / r1
        r2 = diagj( k ) / r2
        a( k-1, k-1 ) = c * r1
        a( k, k-1 ) = s * r1
        a( k-1, k ) = -s * r2
        a( k, k ) = c * r2
     end if

     !     Store details of the interchanges in IPIV
     if ( kstep .eq. 1 ) then
        ipiv( k ) = kp
     else
        ipiv( k ) = -kp
        ipiv( k-1 ) = -imax
     end if

     !     Decrease K and return to the start of the main loop
     k = k - kstep
     go to 10
  else
     !     Factorize A as L*D*L' using the lower triangle of A
     !
     !     K is the main loop index, increasing from 1 to N in steps of 1 or 2
     K = 1
40   continue
     !     If K > N, exit from loop
     if ( k .gt. n ) go to 70

     !     Determine rows and columns to be interchanged and whether
     !     a 1-by-1 or 2-by-2 pivot block will be used

     !     KDIAG is the index of the largest diagonal element, and
     !     DIAMAX is its absolute value
     kdiag = ixamax( n-k+1, a( k, k ), lda + 1 ) + k - 1
     diamax = abs( a( kdiag, kdiag ) )

     !     IMAX and JMAX are the row- and column-indices of the largest
     !     off-diagonal element, and OFFMAX is its absolute value
     offmax = zero
     if ( k .lt. n ) then
        do j = k, n - 1
           do i = j + 1, n
              temp = abs( a( i, j ) )
              if ( temp .gt. offmax ) then
                 offmax = temp
                 imax = i
                 jmax = j
              end if
           end do
        end do
     end if

     if ( max( diamax, offmax ) .eq. zero ) then
        !     The rest of the matrix is zero: set INFO and return
        info = k
        go to 70
     end if

     if ( diamax .ge. alpha * offmax ) then
        !     Use 1-by-1 pivot block
        kstep = 1
        kp = kdiag
     else
        !     Use 2-by-2 pivot block
        kstep = 2
        kp = jmax
     end if

     if ( kp .ne. k ) then
        !     Interchange rows and columns K and KP in the trailing
        !     submatrix A(k:n,k:n)
        if ( kp .lt. n ) call xswap( n-kp, a( kp+1, k ), 1, a( kp+1, kp ), 1 )
        call xswap( kp-k-1, a( k+1, k ), 1, a( kp, k+1 ), lda )
        t = a( k, k )
        a( k, k ) = a( kp, kp )
        a( kp, kp ) = t
     end if

     if ( ( kstep .eq. 2 ) .and. ( imax .ne. k+1 ) ) then
        !     Interchange rows and columns K+1 and IMAX in the trailing
        !     submatrix A(k:n,k:n)
        if ( imax .lt. n ) call xswap( n-imax, a( imax+1, k+1 ), 1, a( imax+1, imax ), 1 )
        call xswap( imax-k-2, a( k+2, k+1 ), 1, a( imax, k+2 ), lda )
        t = a( k+1, k+1 )
        a( k+1, k+1 ) = a( imax, imax )
        a( imax, imax ) = t
        t = a( imax, k )
        a( imax, k ) = a( k+1, k )
        a( k+1, k ) = t
     end if

     !     Update the trailing submatrix
     if ( kstep .eq. 1 ) then
        !     1-by-1 pivot block D(k): column k now holds
        !
        !     W(k+1:n,k) = U(k+1:n,k) * sqrt( abs(D(k)) ),
        !     W(k,k) = U(k,k) * sqrt( abs(D(k)) ) * J(k),
        !
        !     where U(k:n,k) is the k-th column of U, and
        !     J(k) = sign( D(k) )
        !
        !     Perform a rank-1 update of A(k+1:n,k+1:n) as
        !
        !     A := A - U(k)*J(k,k)*U(k)' = A - W(k)*1/D(k)*W(k)'
        r1 = one / a( k, k )
        if ( k .lt. n ) call xsyr( uplo, n-k, -r1, a( k+1, k ), 1, a( k+1, k+1 ), lda )

        !     Compute the k-th diagonal element of the matrix J
        if ( r1 .ge. zero ) then
           diagj( k ) = one
        else
           diagj( k ) = -one
        end if

        !     Store U(k) in column k
        r1 = sqrt( abs( r1 ) )
        !     Singer modif: IF k < n test added!
        if ( k .lt. n ) call xscal( n-k, r1, a( k+1, k ), 1 )
        a( k, k ) = diagj( k ) / r1
     else
        !     2-by-2 pivot block D(k): let
        !
        !     D(k) = Q(k)**T * X(k) * Q(k)
        !
        !     be the eigendecomposition of D(k), X(k) = diag(R1,R2).
        !     Columns k and k-1 now hold
        !
        !     ( W(k+2:n,k) W(k+2:n,k+1) ) =
        !
        !     ( U(k+2:n,k) U(k+2:n,k+1) )*sqrt(abs(X(k)))*Q(k)**T,
        !
        !     W(k:k+1,k:k+1) =
        !
        !     U(k:k+1,k:k+1)*Q(k)*inv(sqrt(abs(X(k))))*J(k),
        !
        !     where U(k) and U(k+1) are the k-th and (k+1)-st columns
        !     of U, and J(k) = diag( sign(R1), sign(R2) ).
        !
        !     Perform a rank-2 update of A(k+2:n,k+2:n) as
        !
        !     A := A - ( U(k) U(k+1) )*J(k)*( U(k) U(k+1) )'
        !     = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )'
        !
        !     Convert this to two rank-1 updates by using the eigen-
        !     decomposition of D(k)
        call xlaev2( a( k, k ), a( k+1, k ), a( k+1, k+1 ), r1, r2, c, s )
        r1 = one / r1
        r2 = one / r2
        call xrot( n-k-1, a( k+2, k ), 1, a( k+2, k+1 ), 1, c, s )
        call xsyr( uplo, n-k-1, -r1, a( k+2, k ), 1, a( k+2, k+2 ), lda )
        call xsyr( uplo, n-k-1, -r2, a( k+2, k+1 ), 1, a( k+2, k+2 ), lda )

        !     Compute the k-th and (k+1)-st diagonal element of the matrix J
        if ( r1 .ge. zero ) then
           diagj( k ) = one
           diagj( k+1 ) = -one
        else
           diagj( k ) = -one
           diagj( k+1 ) = one
        end if

        !     Store U(k) and U(k+1) in columns k and k+1
        r1 = sqrt( abs( r1 ) )
        r2 = sqrt( abs( r2 ) )
        call xscal( n-k-1, r1, a( k+2, k ), 1 )
        call xscal( n-k-1, r2, a( k+2, k+1 ), 1 )

        r1 = diagj( k ) / r1
        r2 = diagj( k+1 ) / r2
        a( k, k ) = c * r1
        a( k, k+1 ) = -s * r2
        a( k+1, k ) = s * r1
        a( k+1, k+1 ) = c * r2
     end if

     !     Store details of the interchanges in IPIV
     if ( kstep .eq. 1 ) then
        ipiv( k ) = kp
     else
        ipiv( k ) = -kp
        ipiv( k+1 ) = -imax
     end if

     !     Increase K and return to the start of the main loop
     k = k + kstep
     go to 40
  end if

70 continue

end subroutine xsyjf2
