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

#ifndef FINT32
typedef int64_t fint;
typedef uint64_t fnat;
#ifndef FINT_C
#define FINT_C(x) INT64_C(x)
#else /* FINT_C */
#error FINT_C already defined
#endif /* ?FINT_C */
#ifndef FNAT_C
#define FNAT_C(x) UINT64_C(x)
#else /* FNAT_C */
#error FNAT_C already defined
#endif /* ?FNAT_C */
#else /* FINT32 */
typedef int32_t fint;
typedef uint32_t fnat;
#ifndef FINT_C
#define FINT_C(x) INT32_C(x)
#else /* FINT_C */
#error FINT_C already defined
#endif /* ?FINT_C */
#ifndef FNAT_C
#define FNAT_C(x) UINT32_C(x)
#else /* FNAT_C */
#error FNAT_C already defined
#endif /* ?FNAT_C */
#endif /* ?FINT32 */

#ifndef F0
#define F0(i, j, ldA) (((j) * (ldA)) + (i))
#endif /* !F0 */

#ifndef F1
#define F1(i, j, ldA) ((((j)-1) * (ldA)) + ((i)-1))
#endif /* !F1 */

#endif /* !COMMON_H */
