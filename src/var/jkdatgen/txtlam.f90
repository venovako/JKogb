subroutine txtlam(n, iseed, lam, npos, info)

  implicit none

  integer, intent(in) :: n
  integer, intent(inout) :: iseed(4)
  real(kind=8), intent(out) :: lam(n)
  integer, intent(out) :: npos, info

  real(kind=8), parameter :: zero = +0.0e+0_8

  integer :: i
  logical :: fex
  character(len=256) :: fn

  external :: seedok

  call seedok(iseed, i)

  if (n .lt. 0) then
    info = -1
  else if (i .ne. 0) then
    info = -2
  else if (command_argument_count() .ne. 7) then
    info = -3
  else
    info = 0
  end if
  if (info .ne. 0) return

  call get_command_argument(7, fn)
  if (len_trim(fn) .eq. 0) then
    info = 1
    return
  end if

  inquire(file=fn, exist=fex)
  if (.not. fex) then
    info = 2
    return
  end if

  open(unit=1,file=fn,action='READ',status='OLD')

  npos = 0
  do i = 1, n
    read (1,*) lam(i)
    if (lam(i) .gt. zero) npos = npos + 1
  end do

  close(unit=1)

end subroutine txtlam
