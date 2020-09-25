#ifndef XTRANSF_H
#define XTRANSF_H
#ifdef USE_EXTENDED
#include "common.h"
extern void zhsvd2_(dcomplex A[static restrict 2][2], const fint J[static restrict 2], dcomplex U[static restrict 2][2], dcomplex Z[static restrict 2][2], fint info[static restrict 1]);
#endif /* USE_EXTENDED */
#endif /* !XTRANSF_H */
