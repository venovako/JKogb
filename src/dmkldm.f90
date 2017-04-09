pure subroutine dmkldm(m, lda)
  implicit none

  integer, intent(in) :: m
  integer, intent(out) :: lda

  integer :: r

  lda = m
  if (m .lt. 0) return

  r = mod(lda, 8)
  if (r .ne. 0) lda = lda + (8 - r)
end subroutine dmkldm
