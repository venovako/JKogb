SUBROUTINE DKOGH2(UPPER, F, G, H, TOL, RHSIDE, CU, SU, CV, SV, F1, H1, INFO)
  IMPLICIT NONE

  LOGICAL, INTENT(IN) :: UPPER
  DOUBLE PRECISION, INTENT(IN) :: F, G, H, TOL

  LOGICAL, INTENT(OUT) :: RHSIDE
  DOUBLE PRECISION, INTENT(OUT) :: CU, SU, CV, SV, F1, H1
  INTEGER, INTENT(OUT) :: INFO

  DOUBLE PRECISION :: TH2, TH, TG, Q

  INFO = 0
  RHSIDE = .TRUE.

  IF (UPPER) THEN
     TH2 = (TWO * F * G) / (F * F + G * G + H * H)
     TH = TH2 / (ONE + SQRT(C_FMA(-TH2, TH2, ONE)))
     TG = C_FMA(F, TH, -G) / H
     CU = ONE / SQRT(C_FMA(TG, TG, ONE))
     SU = CU * TG
     CV = ONE / SQRT(C_FMA(-TH, TH, ONE))
     SV = CV * TH

     Q = CU / CV
     F1 = F * Q
     H1 = H / Q
  ELSE
     TH2 = (TWO * G * H) / (F * F + G * G + H * H)
     TH = TH2 / (ONE + SQRT(C_FMA(-TH2, TH2, ONE)))
     TG = C_FMA(-H, TH, G) / F
     CU = ONE / SQRT(C_FMA(TG, TG, ONE))
     SU = CU * TG
     CV = ONE / SQRT(C_FMA(-TH, TH, ONE))
     SV = CV * TH

     Q = CU / CV
     F1 = F / Q
     H1 = H * Q
  END IF
END SUBROUTINE DKOGH2
