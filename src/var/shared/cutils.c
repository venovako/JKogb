#ifdef USE_INTEL
#include <mathimf.h>
#else /* !USE_INTEL */
#ifdef __cplusplus
#include <cmath>
#else /* !__cplusplus */
#include <math.h>
#endif /* ?__cplusplus */
#endif /* ?USE_INTEL */

#include "cutils.h"

static inline long double XASUM2(const double a, const double b) {
  const double a_ = fabs(a);
  const double b_ = fabs(b);
  const long double x = (long double)a_;
  const long double y = (long double)b_;
  return ((b_ <= a_) ? fmal(x, x, (y * y)) : fmal(y, y, (x * x)));
}

static inline long double XASUM4(const double a, const double b, const double c, const double d)
{
  long double z[4];
  double x[4];
  double y[4];

  x[0] = fabs(a);
  x[1] = fabs(b);
  x[2] = fabs(c);
  x[3] = fabs(d);

  if (x[0] <= x[1]) {
    y[0] = x[0];
    y[1] = x[1];
  }
  else {
    y[0] = x[1];
    y[1] = x[0];
  }

  if (x[2] <= x[3]) {
    y[2] = x[2];
    y[3] = x[3];
  }
  else {
    y[2] = x[3];
    y[3] = x[2];
  }

  if (y[0] <= y[2]) {
    x[0] = y[0];
    x[1] = y[2];
  }
  else {
    x[0] = y[2];
    x[1] = y[0];
  }

  if (y[1] <= y[3]) {
    x[2] = y[1];
    x[3] = y[3];
  }
  else {
    x[2] = y[3];
    x[3] = y[1];
  }

  z[0] = (long double)(x[0]);
  if (x[1] <= x[2]) {
    z[1] = (long double)(x[1]);
    z[2] = (long double)(x[2]);
  }
  else {
    z[1] = (long double)(x[2]);
    z[2] = (long double)(x[1]);
  }
  z[3] = (long double)(x[3]);

  return fmal(z[3], z[3], fmal(z[2], z[2], fmal(z[1], z[1], (z[0] * z[0]))));
}

double HYPOTwX87(const double a, const double b)
{
  return (double)sqrtl(XASUM2(a, b));
}

double DASUM2(const double a, const double b)
{
  return (double)XASUM2(a, b);
}

double DASUM4(const double a, const double b, const double c, const double d)
{
  return (double)XASUM4(a, b, c, d);
}
