#ifndef DJK2_H
#define DJK2_H

#include "common.h"

#ifdef USE_EXTENDED
extern void dhsvd2_(double A[static 2][2], const fint J[static 2], double U[static 2][2], double Z[static 2][2], fint info[static 1]);
#endif /* USE_EXTENDED */

#endif /* !DJK2_H */
