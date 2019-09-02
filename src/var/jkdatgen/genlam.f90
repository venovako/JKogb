subroutine genlam(n, iseed, lam, npos, info)

  implicit none

  integer, intent(in) :: n
  integer, intent(inout) :: iseed(4)
  real(kind=8), intent(out) :: lam(n)
  integer, intent(out) :: npos, info

  real(kind=8), parameter :: zero = +0.0e+0_8, one = +1.0e+0_8

  integer :: idist, i
  real(kind=8) :: eps, scal
  character(len=256) :: cas

  real(kind=8), external :: dlarnd
  external :: seedok

  call seedok(iseed, i)

  if (n .lt. 0) then
    info = -1
  else if (i .ne. 0) then
    info = -2
  else if (command_argument_count() .ne. 9) then
    info = -3
  else
    info = 0
  end if
  if (info .ne. 0) return

  call get_command_argument(7, cas)
  read (cas,*) idist
  if ((idist .lt. 1) .or. (idist .gt. 3)) then
    info = 1
    return
  end if

  call get_command_argument(8, cas)
  read (cas,*) eps
  if (eps .le. zero) then
    info = 2
    return
  end if

  call get_command_argument(9, cas)
  read (cas,*) scal
  if (scal .eq. zero) then
    info = 3
    return
  end if

  npos = 0
  i = 1

  do while (i .le. n)
    lam(i) = dlarnd(idist, iseed)
    if (abs(lam(i)) .gt. eps) then
      if (lam(i) .gt. zero) npos = npos + 1
      i = i + 1
    end if
  end do

  if (scal .ne. one) then
    do i = 1, n
      lam(i) = lam(i) * scal
    end do
  end if

end subroutine genlam
