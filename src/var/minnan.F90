program minnan
  double precision :: x
  x = 0d0
  x = x / x ! x is NaN
  print *, 'MIN(x,1)=', min(x,1d0)
  print *, 'MIN(1,x)=', min(1d0,x)
end program minnan
