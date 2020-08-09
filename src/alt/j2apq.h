#ifndef J2APQ_H
#define J2APQ_H
#include "wpqb.h"

typedef struct {
  uint32_t n_a, n_t;
  wpqb *a;
  uint16_t *p, *q;
} apq;

extern int apq_init(apq o[static 1], const uint16_t n, const fint *const j);
extern void apq_free(apq o[static 1]);
#endif /* !J2APQ_H */
