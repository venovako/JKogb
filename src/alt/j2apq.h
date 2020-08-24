#ifndef J2APQ_H
#define J2APQ_H
#include "wpqb.h"

typedef struct {
  uint32_t n_a, n_t;
  wpqb *a;
  uint16_t *p, *q;
} apq;

extern apq *apq_create(const uint16_t n, const fint *const j);
extern apq *apq_free(apq *const o);
#endif /* !J2APQ_H */
