PURE SUBROUTINE CA(B, N, X, Y, LDA)
  IMPLICIT NONE
  COMPLEX(KIND=DWP), INTENT(IN) :: B(2,2)
  INTEGER, INTENT(IN) :: N, LDA
  COMPLEX(KIND=DWP), INTENT(INOUT) :: X(*), Y(*)

  COMPLEX(KIND=DWP) :: XX, YY
  INTEGER :: I, J

  I = 1

  ! DO J = 1, N
  !    XX = X(I) + B(1,2) * Y(I)
  !    YY = B(2,1) * X(I) + Y(I)
  !    X(I) = B(1,1) * XX
  !    Y(I) = B(2,2) * YY
  !    I = I + LDA
  ! END DO

  DO J = 1, N
     XX = ZFMA(B(1,2), Y(I), X(I))
     YY = ZFMA(B(2,1), X(I), Y(I))
     X(I) = XX / REAL(B(1,1))
     Y(I) = YY / REAL(B(2,2))
     I = I + LDA
  END DO
END SUBROUTINE CA
