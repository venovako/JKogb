#ifndef J2PQ_H
#define J2PQ_H
#include "common.h"

typedef struct {
  uint32_t n_a, n_t;
  uint16_t *p, *q;
} pq;

extern void pq_init(pq o[static 1], const uint16_t n, const int j[static 1]);
extern void pq_free(pq o[static 1]);
#endif /* !J2PQ_H */
