subroutine gendat(lgix, seed, n, lam, iseed, nrank, a, g, ldg, ipiv, jvec, nplus, ipl, invp, info)

  implicit none

  integer, intent(in) :: lgix, seed, n, ldg
  integer, intent(out) :: iseed(4), nrank, ipiv(n), jvec(n), nplus, ipl(n), invp(n), info(2)
  real(kind=8), intent(out) :: lam(n), a(ldg,n), g(ldg,n)

  integer :: jseed(4), npos

  external :: seedix, txtlam, genlam, xa2g

  info(1) = 0
  info(2) = 0

  call seedix(seed, iseed, info(2))
  if (info(2) .ne. 0) then
     info(1) = 10
     return
  end if
  jseed = iseed

  select case (lgix)
  case (1)
    call txtlam(n, jseed, lam, npos, info(2))
    if (info(2) .ne. 0) info(1) = 21
  case (2)
    call genlam(n, jseed, lam, npos, info(2))
    if (info(2) .ne. 0) info(1) = 22
  case default
    info(1) = 20
    info(2) = lgix
  end select
  if (info(1) .ne. 0) return

  call xa2g(n, lam, jseed, nrank, a, g, ldg, ipiv, jvec, nplus, ipl, invp, info)
  if ((info(1) .eq. 0) .and. (nplus .ne. npos)) then
     info(1) = 30
     info(2) = npos
  end if

end subroutine gendat
