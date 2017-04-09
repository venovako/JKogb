subroutine dkogh2(upper, f, g, h, tol, rhside, cu, su, cv, sv, f1, h1, info)
  implicit none

  logical, intent(in) :: upper
  double precision, intent(in) :: f, g, h, tol

  logical, intent(out) :: rhside
  double precision, intent(out) :: cu, su, cv, sv, f1, h1
  integer, intent(out) :: info

  double precision :: th2, th, tg, q
  logical :: flag

  info = 0
  rhside = .true.

  if (upper) then
     th2 = (TWO * f * g) / (f * f + g * g + h * h)
     th = th2 / (ONE + sqrt(c_fma(-th2, th2, ONE)))
     tg = c_fma(f, th, -g) / h
     cu = ONE / sqrt(c_fma(tg, tg, ONE))
     su = cu * tg
     cv = ONE / sqrt(c_fma(-th, th, ONE))
     sv = cv * th

     q = cu / cv
     f1 = f * q
     h1 = h / q
  else
     th2 = (TWO * g * h) / (f * f + g * g + h * h)
     th = th2 / (ONE + sqrt(c_fma(-th2, th2, ONE)))
     tg = c_fma(-h, th, g) / f
     cu = ONE / sqrt(c_fma(tg, tg, ONE))
     su = cu * tg
     cv = ONE / sqrt(c_fma(-th, th, ONE))
     sv = cv * th

     q = cu / cv
     f1 = f / q
     h1 = h * q
  end if
end subroutine dkogh2
