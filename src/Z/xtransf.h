#ifndef XTRANSF_H
#define XTRANSF_H
#ifdef USE_EXTENDED
#include "common.h"
extern void zhsvd2_(dcomplex A[static 2][2], const fint J[static 2], dcomplex U[static 2][2], dcomplex Z[static 2][2], fint info[static 1]);
#endif /* USE_EXTENDED */
#endif /* !XTRANSF_H */
