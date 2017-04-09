  interface
     real(c_double) function c_fma(x, y, z) bind(c,name='fma')
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), value :: x, y, z
     end function c_fma
  end interface
