program itrue
  implicit none
  logical :: l
  integer :: i
  l = .true.
  i = transfer(l,i)
  print *, i
end program itrue
