SUBROUTINE DKOGT2(UPPER, F, G, H, TOL, RHSIDE, CU, SU, CV, SV, F1, H1, INFO)
  IMPLICIT NONE

  LOGICAL, INTENT(IN) :: UPPER
  DOUBLE PRECISION, INTENT(IN) :: F, G, H, TOL

  LOGICAL, INTENT(OUT) :: RHSIDE
  DOUBLE PRECISION, INTENT(OUT) :: CU, SU, CV, SV, F1, H1
  INTEGER, INTENT(OUT) :: INFO

  INFO = 0

  IF (UPPER) THEN
     IF (F .GE. H) THEN
        RHSIDE = .FALSE.
        CALL DKOGUL(F, G, H, TOL, CU, SU, CV, SV, F1, H1)
     ELSE
        RHSIDE = .TRUE.
        CALL DKOGUL(H, -G, F, TOL, CV, SV, CU, SU, H1, F1)
     END IF
  ELSE
     IF (F .GE. H) THEN
        RHSIDE = .TRUE.
        CALL DKOGUL(F, G, H, TOL, CV, SV, CU, SU, F1, H1)
     ELSE
        RHSIDE = .FALSE.
        CALL DKOGUL(H, -G, F, TOL, CU, SU, CV, SV, H1, F1)
     END IF
  END IF
END SUBROUTINE DKOGT2