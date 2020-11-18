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

double HYPOTwX87(const double a, const double b)
{
  long double x = (long double)a;
  long double y = (long double)b;
  x = x * x;
  y = y * y;
  return (double)sqrtl(x + y);
}
