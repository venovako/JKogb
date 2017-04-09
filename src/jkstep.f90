pure logical function ldorot(f, g, h, tol, hyp)
  implicit none

  double precision, intent(in) :: f, g, h, tol
  logical, intent(in) :: hyp

  double precision :: ga

  ga = abs(g)
  ldorot = (f .ne. (f + ga)) .or. (h .ne. (h + ga))
end function ldorot

subroutine jkstep(upper, i, j, m, n, G, ldg, U, ldu, V, ldv, jvec, nplus, tol, sweep, swhalf, step, hyp, info)
  implicit none

  logical, intent(in) :: upper
  integer, intent(in) :: i, j, m, n, ldg, ldu, ldv, jvec(n), nplus, sweep, swhalf
  double precision, intent(in) :: tol
  double precision, intent(inout) :: G(ldg, n), U(ldu, m), V(ldv, n)
  integer, intent(inout) :: step
  logical, intent(out) :: hyp
  integer, intent(out) :: info(2)

  double precision :: f0, g0, h0, cu, su, cv, sv, f1, h1, paramu(5), paramv(5)
  logical :: rhside
  integer :: instat

  external :: drotm

  info = 0
  step = step + 1
  hyp = jvec(i) .ne. jvec(j)

  f0 = G(i, i)
  h0 = G(j, j)

  if (upper) then
     g0 = G(i, j)
  else
     g0 = G(j, i)
  end if

  if (g0 .eq. ZERO) return

  if (ldorot(f0, g0, h0, tol, hyp)) then
     info(2) = 1
     paramu(1) = MONE
     paramv(1) = MONE

     if (hyp) then
        call dkogh2(upper, f0, g0, h0, tol, rhside, cu, su, cv, sv, f1, h1, instat)
     else
        call dkogt2(upper, f0, g0, h0, tol, rhside, cu, su, cv, sv, f1, h1, instat)
     end if

     paramu(2) =  cu
     paramu(3) = -su
     paramu(4) =  su
     paramu(5) =  cu

     if (hyp) then
        paramv(2) =  cv
        paramv(3) = -sv
        paramv(4) = -sv
        paramv(5) =  cv
     else
        paramv(2) =  cv
        paramv(3) = -sv
        paramv(4) =  sv
        paramv(5) =  cv
     end if

     if (rhside) then
        call drotm(m, G(1, i), 1, G(1, j), 1, paramv)
        call drotm(n, G(i, 1), ldg, G(j, 1), ldg, paramu)
     else
        call drotm(n, G(i, 1), ldg, G(j, 1), ldg, paramu)
        call drotm(m, G(1, i), 1, G(1, j), 1, paramv)
     end if

     call drotm(m, U(i, 1), ldu, U(j, 1), ldu, paramu)
     call drotm(n, V(1, i), 1, V(1, j), 1, paramv)

     G(i, i) = f1
     G(j, j) = h1
  else
     info(2) = -1
  end if
  G(i, j) = ZERO
  G(j, i) = ZERO
end subroutine jkstep
