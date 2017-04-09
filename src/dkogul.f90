subroutine dkogul(f, g, h, tol, c1, s1, c2, s2, f1, h1)
  implicit none

  double precision, intent(in) :: f, g, h, tol
  double precision, intent(out) :: c1, s1, c2, s2, f1, h1

  double precision :: d, e, zeta, mu, mu2, rho, rho2, alpha, t1, t2, x, q

  if (abs(g) .gt. f) then
     !d = g + ((f + h) / g) * (f - h)
     d = c_fma((f + h) / g, f - h, g)
     e = TWO * h
  else
     d = (f - h) / g + g / (f + h)
     e = TWO * h / (f + h)
  end if

  if (abs(d) .le. e) then
     zeta = d / e
     !mu = ONE + zeta * zeta
     mu = c_fma(zeta, zeta, ONE)
     rho = sqrt(mu)
     alpha = ONE + abs(zeta) / rho
     c1 = sqrt(HALF * alpha)
     s1 = sign(ONE / sqrt(TWO * alpha * mu), zeta)
     t1 = sign(ONE / (abs(zeta) + rho), zeta)
  else
     if (abs(d) * tol .gt. e) then
        c1 = ONE
        t1 = (HALF * e) / d
        s1 = t1
        zeta = ZERO
        alpha = TWO
     else
        zeta = e / d
        !mu = ONE + zeta * zeta
        mu = c_fma(zeta, zeta, ONE)
        rho = sqrt(mu)
        alpha = ONE + ONE / rho
        c1 = sqrt(HALF * alpha)
        s1 = zeta / sqrt(TWO * alpha * mu)
        t1 = zeta / (ONE + rho)
     end if
  end if

  !x = h * t1 + g
  x = c_fma(h, t1, g)
  if (abs(x) .le. f) then
     t2 = x / f
     !mu2 = ONE + t2 * t2
     mu2 = c_fma(t2, t2, ONE)
     rho2 = sqrt(mu2)
     c2 = ONE / rho2
     s2 = t2 / rho2
     q = sqrt(HALF * alpha * mu2)
  else
     t2 = f / x
     !mu2 = ONE + t2 * t2
     mu2 = c_fma(t2, t2, ONE)
     rho2 = sqrt(mu2)
     s2 = sign(ONE / rho2, t2)
     c2 = abs(t2) / rho2
     q = sqrt(HALF * alpha * mu2) / abs(t2)
  end if

  f1 = f * q
  h1 = h / q
end subroutine dkogul
