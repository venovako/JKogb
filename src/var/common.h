#ifndef COMMON_H
#define COMMON_H

#ifdef __ICC
#include <mathimf.h>
#else /* !__ICC */
#include <complex.h>
#include <math.h>
#endif /* __ICC */

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <fenv.h>
#include <float.h>
#include <limits.h>
#include <signal.h>
#include <stdalign.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>
#include <time.h>

/* Intel 80-bit extended floating-point value stored in the lowest 10 bytes of a 16-byte variable */
typedef long double extended;

#endif /* !COMMON_H */
