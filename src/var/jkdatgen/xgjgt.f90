subroutine xgjgt( info, n, nrank, g, lda, ipiv, jvec, infog )

  !     Purpose
  !     =======
  !
  !     GJGT computes the indefinite factorization from XSYTRF.
  !     Returns G in the LOWER triangle of G.

  implicit none

  integer, intent(in) :: info
  integer, intent(in) :: n
  integer, intent(inout) :: nrank
  real(kind=16), dimension(lda, n), intent(inout) :: g
  integer, intent(in) :: lda
  integer, dimension(n), intent(inout) :: ipiv
  integer, dimension(n), intent(out) :: jvec
  integer, intent(out) :: infog

  real(kind=16), parameter :: zero = +0.0e+0_16

  integer :: nr2, k, kp
  real(kind=16) :: aa, bb, cc, cs1, sn1, d, e

  external :: xlaev2, xlaset, xrot, xscal, xswap

  infog = 0
  if ( info .lt. 0 ) return

  if ( nrank .lt. 0 ) then
     infog = -1
  else if ( n .lt. 0 ) then
     infog = -2
  else if ( lda .lt. max( 1, n ) ) then
     infog = -4
  end if

  nrank = n
  if ( info .gt. 0 ) nrank = info - 1
  if ( nrank .eq. 0 ) return

  k = 1
10 continue
  if ( k .gt. nrank ) go to 20
  if ( ipiv( k ) .gt. 0 ) then
     d = g( k, k )
     if ( d .lt. zero ) then
        jvec( k ) = -1
     else
        jvec( k ) = 1
     end if
     d = sqrt( abs( d ) )
     if ( k .lt. n ) then
        call xscal( n-k, d, g( k+1, k ), 1 )
     end if
     g( k, k ) = d
     k = k + 1
  else if ( ipiv( k ) .lt. 0 ) then
     aa = g( k, k )
     bb = g( k+1, k )
     cc = g( k+1, k+1 )
     call xlaev2( aa, bb, cc, d, e, cs1, sn1 )
     if ( d .lt. zero ) then
        jvec( k ) = -1
        jvec( k+1 ) = 1
     else
        jvec( k ) = 1
        jvec( k+1 ) = -1
     end if
     d = sqrt( abs( d ) )
     e = sqrt( abs( e ) )
     if ( k+1 .lt. n ) then
        call xrot( n-k-1, g( k+2, k ), 1, g( k+2, k+1 ), 1, cs1, sn1 )
        call xscal( n-k-1, d, g( k+2, k ), 1 )
        call xscal( n-k-1, e, g( k+2, k+1 ), 1 )
     end if
     g( k, k ) = d * cs1
     g( k+1, k+1 ) = e * cs1
     g( k, k+1 ) = -e * sn1
     g( k+1, k ) = d * sn1
     k = k + 2
  end if
  go to 10

20 continue
  !     Make P*H*P' = U*J*U'   or  P*H*P' = L*J*L'

  !     Create the factor G such that  P*A*P' = G*J*G', where G is
  !     lower block triangular matrix, and J is as above.

  !     Empty the the the array G above the second superdiagonal
  nr2 = nrank - 2
  if ( nr2 .gt. 0 ) call xlaset( 'U', nr2, nr2, zero, zero, g( 1, 3 ), lda )

  !     Set the element G( 1, 2 ) or G( 2, 3 ) to zero if necessary
  if ( ipiv( 1 ) .gt. 0 ) then
     if ( nrank .gt. 1 ) g( 1, 2 ) = zero
     !     K is the loop index, increasing from 2 or 3 to NRANK in steps
     !     of 1 or 2, depending on the size of the diagonal blocks.
     k = 2
  else
     if ( nrank .gt. 2 ) g( 2, 3 ) = zero
     k = 3
  end if

30 continue
  !     If K > NRANK, exit from loop.
  if ( k .gt. nrank ) go to 40

  if ( ipiv( k ) .gt. 0 ) then
     !     1 x 1 diagonal block

     !     Interchange rows K and IPIV(K).
     kp = ipiv( k )
     if ( kp .ne. k ) call xswap( k-1, g( k, 1 ), lda, g( kp, 1 ), lda )

     !     Set the element G( k, k+1 ) of the upper triangle to zero
     if ( nrank .gt. k ) g( k, k+1 ) = zero
     k = k + 1
  else
     !     2 x 2 diagonal block

     !     NO Interchange for row K.
     !     Interchange rows K+1 and -IPIV(K+1).
     kp = -ipiv( k+1 )
     if ( kp .ne. k+1 ) call xswap( k-1, g( k+1, 1 ), lda, g( kp, 1 ), lda )

     !     Set the element G( k+1, k+2 ) of the upper triangle to zero
     if ( nrank .gt. k+1 ) g( k+1, k+2 ) = zero
     k = k + 2
  end if
  go to 30

40 continue

end subroutine xgjgt
