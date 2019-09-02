subroutine jpart(nrow, ncolr, g, ldg, jvec, nplus, ipl, invp)

  !     Purpose
  !     =======
  !
  !     Transforms the G*J*G^T factorization into G_1*J_part*G_1^T
  !     factorization, with J_part partitioned as  J_part = ( I, -I ).
  !     Reorders the columns of G and the elements of JVEC.
  !     NPLUS is the number of elements in JVEC equal to 1.

  implicit none

  integer, intent(in) :: nrow, ncolr, ldg
  real(kind=16), dimension(ldg, nrow), intent(inout) :: g
  integer, dimension(nrow), intent(inout) :: jvec
  integer, intent(out) :: nplus, ipl(nrow), invp(nrow)

  integer :: i, iplus, iminus, ip, jtemp

  external :: xswap

  !     Count columns with JVEC( I ) = 1.
  nplus = 0
  do i = 1, ncolr
     if ( jvec( i ) .eq. 1 ) nplus = nplus + 1
  end do

  !     Set permutation IPL, where IPL( I ) holds the current place
  !     of the final I-th column.
  !     The following algorithm preserves the relative order of columns
  !     with the same sign in JVEC.
  iplus = 0
  iminus = nplus
  do i = 1, ncolr
     if ( jvec( i ) .eq. 1 ) then
        iplus = iplus + 1
        ipl( iplus ) = i
     else
        iminus = iminus + 1
        ipl( iminus ) = i
     end if
  end do
  do i = ncolr + 1, nrow
     ipl( i ) = i
  end do

  !     Invert the permutation IPL and store it in INVP.
  do i = 1, nrow
     invp( ipl( i ) ) = i
  end do

  !     Early return - all JVEC( I ) have the same sign.
  if ( ( nplus .eq. 0 ) .or. ( nplus .eq. ncolr ) ) goto 1

  do i = 1, ncolr
     !     Swap columns G( I ) and G( IPL( I ) ).
     !     Also swap the corresponding elements in JVEC.
     if ( ipl( i ) .ne. i ) then
        ip = ipl( i )

        call xswap( nrow, g( 1, i ), 1, g( 1, ip ), 1 )
        jtemp = jvec( i )
        jvec( i ) = jvec( ip )
        jvec( ip ) = jtemp

        invp( ip ) = invp( i )
        ipl( invp( i ) ) = ip
     end if
     !     Not necessary to set:
     !     IPL( I ) = I
     !     INVP( I ) = I
  end do

1 do i = ncolr + 1, nrow
     jvec( i ) = 0
  end do

end subroutine jpart
