#ifndef COMMON_H
#define COMMON_H

#if (defined(__ICC) || defined(__INTEL_CLANG_COMPILER) || defined(__INTEL_LLVM_COMPILER))
#include <mathimf.h>
#else /* !__ICC */
#include <complex.h>
#include <math.h>
#endif /* ?__ICC */

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
#ifndef USE_NVIDIA
#include <tgmath.h>
#endif /* !USE_NVIDIA */
#include <time.h>

// Intel 80-bit extended floating-point value stored in the lowest 10 bytes of a 16-byte variable
typedef long double extended;
typedef double complex dcomplex;
typedef long double complex xcomplex;

#ifndef CMPLXF
#define CMPLXF(r,i) ((float)(r) + I * (float)(i))
#endif /* !CMPLXF */
#ifndef CMPLX
#define CMPLX(r,i) ((double)(r) + I * (double)(i))
#endif /* !CMPLX */
#ifndef CMPLXL
#define CMPLXL(r,i) ((extended)(r) + I * (extended)(i))
#endif /* !CMPLXL */

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

// sqrtl(LDBL_MAX)
#ifndef TWOF
#define TWOF 1.090748135619415929404E+2466L
#else /* TWOF */
#error TWOF already defined
#endif /* ?TWOF */

#endif /* !COMMON_H */
