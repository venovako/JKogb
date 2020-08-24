#ifndef COMMON_H
#define COMMON_H

#ifdef __ICC
#include <mathimf.h>
#else /* !__ICC */
#ifdef _OPENMP
#error OpenMP variant requires a recent Intel C compiler
#endif /* _OPENMP */
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
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef CMPLXF
#define CMPLXF(r,i) ((float)(r) + I * (float)(i))
#endif /* !CMPLXF */
#ifndef CMPLX
#define CMPLX(r,i) ((double)(r) + I * (double)(i))
#endif /* !CMPLX */
#ifndef CMPLXL
#define CMPLXL(r,i) ((long double)(r) + I * (long double)(i))
#endif /* !CMPLXL */

#ifdef FINT64
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
#else /* !FINT64 */
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
#endif /* ?FINT64 */

// sqrtl(LDBL_MAX)
#ifndef TWOF
#define TWOF 1.090748135619415929404E+2466L
#else /* TWOF */
#error TWOF already defined
#endif /* ?TWOF */

#endif /* !COMMON_H */
