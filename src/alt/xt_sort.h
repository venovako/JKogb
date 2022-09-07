#ifndef XT_SORT_H
#define XT_SORT_H

#include <mkl.h>

#ifdef __cplusplus
extern "C" void xt_sort_(const MKL_INT *const n, long double *const x) throw();
#else /* !__cplusplus */
extern void xt_sort_(const MKL_INT n[static restrict 1], long double x[static restrict 1]);
#endif /* ?__cplusplus */

#endif /* !XT_SORT_H */
