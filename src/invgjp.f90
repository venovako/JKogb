subroutine invgjp(nrow, ncolr, G, ldg, jvec, ipl, invp)
  implicit none

  integer, intent(in) :: nrow, ncolr, ldg
  double precision, intent(inout) :: G(ldg, ncolr)
  integer, intent(inout) :: jvec(ncolr), ipl(ncolr), invp(ncolr)

  integer :: i, ip, jtemp

  external :: dswap

  do i = 1, ncolr
    if (invp(i) .ne. i) then
       ip = invp(i)

       call dswap(nrow, G(1, i), 1, G(1, ip), 1)
       jtemp = jvec(i)
       jvec(i) = jvec(ip)
       jvec(ip) = jtemp
       
       ipl(ip) = ipl(i)
       invp(ipl(i)) = ip
    end if
 end do
end subroutine invgjp
