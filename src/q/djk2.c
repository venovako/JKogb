#include "djk2.h"

int Ut2(double *const U, const size_t ldU)
{
  if (U) {
    if (ldU >= 2u) {
      const double u10 = U[F0(1, 0, ldU)];
      U[F0(1, 0, ldU)] = U[F0(0, 1, ldU)];
      U[F0(0, 1, ldU)] = u10;
      return 0;
    }
    return -2;
  }
  return -1;
}

int JZtJ2(double *const Z, const size_t ldZ, const int *const J)
{
  if (Z) {
    if (ldZ >= 2u) {
      if (J) {
        const int t = Ut2(Z, ldZ);
        if (t)
          return -t;
        switch (J[0]) {
        case 1:
          Z[F0(1, 0, ldZ)] = -Z[F0(1, 0, ldZ)];
          Z[F0(0, 1, ldZ)] = -Z[F0(0, 1, ldZ)];
          break;
        case -1:
          break;
        default:
          return -3;
        }
        switch (J[1]) {
        case 1:
          break;
        case -1:
          Z[F0(1, 0, ldZ)] = -Z[F0(1, 0, ldZ)];
          Z[F0(0, 1, ldZ)] = -Z[F0(0, 1, ldZ)];
          break;
        default:
          return -4;
        }
        return 0;
      }
      return Ut2(Z, ldZ);
    }
    return -2;
  }
  return -1;
}
