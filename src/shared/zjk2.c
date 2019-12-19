#include "zjk2.h"

void zhsvd2_(dcomplex A[static 2][2], const fint J[static 2], dcomplex U[static 2][2], dcomplex Z[static 2][2], fint info[static 1])
{
  *info = FINT_C(0);

  switch (J[0]) {
  case FINT_C(-1):
  case FINT_C(1):
    break;
  default:
    *info = FINT_C(-5);
    return;
  }

  switch (J[1]) {
  case FINT_C(-1):
  case FINT_C(1):
    break;
  default:
    *info = FINT_C(-6);
    return;
  }

  xcomplex A_[2][2] = { { A[0][0], A[0][1] }, { A[1][0], A[1][1] } };
  xcomplex U_[2][2] = { { CMPLXL(1.0L, 0.0L), CMPLXL(0.0L, 0.0L) }, { CMPLXL(0.0L, 0.0L), CMPLXL(1.0L, 0.0L) } };
  xcomplex Z_[2][2] = { { CMPLXL(1.0L, 0.0L), CMPLXL(0.0L, 0.0L) }, { CMPLXL(0.0L, 0.0L), CMPLXL(1.0L, 0.0L) } };
}
