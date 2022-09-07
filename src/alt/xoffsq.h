#ifndef XOFFSQ_H
#define XOFFSQ_H

#include <mkl.h>

#ifdef __cplusplus
extern "C" void xoffsq_(const MKL_INT *const n, const double *const A, const MKL_INT *const ldA, long double *const x, MKL_INT *const info) throw();
#else /* !__cplusplus */
extern void xoffsq_(const MKL_INT n[static 1], const double A[static restrict 1], const MKL_INT ldA[static 1], long double x[static restrict 1], MKL_INT info[static restrict 1]);
#endif /* ?__cplusplus */

#endif /* !XOFFSQ_H */
