  INTERFACE
     FUNCTION C_FMA(X, Y, Z) BIND(C,NAME='fma')
       USE, INTRINSIC :: ISO_C_BINDING
       IMPLICIT NONE
       REAL(c_double), VALUE :: X, Y, Z
       REAL(c_double) :: C_FMA
     END FUNCTION C_FMA
  END INTERFACE
