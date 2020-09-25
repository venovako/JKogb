PURE SUBROUTINE AC(B, M, X, Y)
  IMPLICIT NONE
  COMPLEX(KIND=DWP), INTENT(IN) :: B(2,2)
  INTEGER, INTENT(IN) :: M
  COMPLEX(KIND=DWP), INTENT(INOUT) :: X(M), Y(M)

  COMPLEX(KIND=DWP) :: XX, YY
  INTEGER :: I

  ! DO I = 1, M
  !    XX = X(I) + Y(I) * B(2,1)
  !    YY = X(I) * B(1,2) + Y(I)
  !    X(I) = XX * REAL(B(1,1))
  !    Y(I) = YY * REAL(B(2,2))
  ! END DO

  DO I = 1, M
     XX = ZFMA(Y(I), B(2,1), X(I))
     YY = ZFMA(X(I), B(1,2), Y(I))
     X(I) = XX / REAL(B(1,1))
     Y(I) = YY / REAL(B(2,2))
  END DO
END SUBROUTINE AC
