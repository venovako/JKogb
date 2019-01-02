SUBROUTINE DKOGUL(F, G, H, TOL, C1, S1, C2, S2, F1, H1)
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: F, G, H, TOL
  DOUBLE PRECISION, INTENT(OUT) :: C1, S1, C2, S2, F1, H1

  DOUBLE PRECISION :: D, E, ZETA, MU, MU2, RHO, RHO2, ALPHA, T1, T2, X, Q

  IF (ABS(G) .GT. F) THEN
     ! D = G + ((F + H) / G) * (F - H)
     D = C_FMA((F + H) / G, F - H, G)
     E = TWO * H
  ELSE
     D = (F - H) / G + G / (F + H)
     E = TWO * H / (F + H)
  END IF

  IF (ABS(D) .LE. E) THEN
     ZETA = D / E
     ! MU = ONE + ZETA * ZETA
     MU = C_FMA(ZETA, ZETA, ONE)
     RHO = SQRT(MU)
     ALPHA = ONE + ABS(ZETA) / RHO
     C1 = SQRT(HALF * ALPHA)
     S1 = SIGN(ONE / SQRT(TWO * ALPHA * MU), ZETA)
     T1 = SIGN(ONE / (ABS(ZETA) + RHO), ZETA)
  ELSE
     IF (ABS(D) * TOL .GT. E) THEN
        C1 = ONE
        T1 = (HALF * E) / D
        S1 = T1
        ZETA = ZERO
        ALPHA = TWO
     ELSE
        ZETA = E / D
        ! MU = ONE + ZETA * ZETA
        MU = C_FMA(ZETA, ZETA, ONE)
        RHO = SQRT(MU)
        ALPHA = ONE + ONE / RHO
        C1 = SQRT(HALF * ALPHA)
        S1 = ZETA / SQRT(TWO * ALPHA * MU)
        T1 = ZETA / (ONE + RHO)
     END IF
  END IF

  ! X = H * T1 + G
  X = C_FMA(H, T1, G)
  IF (ABS(X) .LE. F) THEN
     T2 = X / F
     ! MU2 = ONE + T2 * T2
     MU2 = C_FMA(T2, T2, ONE)
     RHO2 = SQRT(MU2)
     C2 = ONE / RHO2
     S2 = T2 / RHO2
     Q = SQRT(HALF * ALPHA * MU2)
  ELSE
     T2 = F / X
     ! MU2 = ONE + T2 * T2
     MU2 = C_FMA(T2, T2, ONE)
     RHO2 = SQRT(MU2)
     S2 = SIGN(ONE / RHO2, T2)
     C2 = ABS(T2) / RHO2
     Q = SQRT(HALF * ALPHA * MU2) / ABS(T2)
  END IF

  F1 = F * Q
  H1 = H / Q
END SUBROUTINE DKOGUL
