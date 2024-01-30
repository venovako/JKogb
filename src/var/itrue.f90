program itrue
  implicit none
  logical :: l
  integer :: i
  l = .false.
  i = transfer(l,i)
  print *, '.false.=', i
  l = .true.
  i = transfer(l,i)
  print *, '.true. =', i
end program itrue
