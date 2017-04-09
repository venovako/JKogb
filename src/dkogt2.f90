subroutine dkogt2(upper, f, g, h, tol, rhside, cu, su, cv, sv, f1, h1, info)
  implicit none

  logical, intent(in) :: upper
  double precision, intent(in) :: f, g, h, tol

  logical, intent(out) :: rhside
  double precision, intent(out) :: cu, su, cv, sv, f1, h1
  integer, intent(out) :: info

  info = 0

  if (upper) then
     if (f .ge. h) then
        rhside = .false.
        call dkogul(f, g, h, tol, cu, su, cv, sv, f1, h1)
     else
        rhside = .true.
        call dkogul(h, -g, f, tol, cv, sv, cu, su, h1, f1)
     end if
  else
     if (f .ge. h) then
        rhside = .true.
        call dkogul(f, g, h, tol, cv, sv, cu, su, f1, h1)
     else
        rhside = .false.
        call dkogul(h, -g, f, tol, cu, su, cv, sv, h1, f1)
     end if
  end if
end subroutine dkogt2
