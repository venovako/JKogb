#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main()
{
  const long double z = sqrtl(LDBL_MAX);
  (void)printf("sqrtl(LDBL_MAX)=%#.21LE, fmal(sqrtl(LDBL_MAX),sqrtl(LDBL_MAX),1)=%#.21LE\n", z, fmal(z, z, 1.0));
  return EXIT_SUCCESS;
}
