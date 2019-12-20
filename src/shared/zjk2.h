#ifndef ZJK2_H
#define ZJK2_H

#include "common.h"

#ifdef USE_EXTENDED
extern void zhsvd2_(dcomplex A[static 2][2], const fint J[static 2], dcomplex U[static 2][2], dcomplex Z[static 2][2], fint info[static 1]);
#endif /* USE_EXTENDED */

#endif /* !ZJK2_H */
