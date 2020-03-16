#ifndef WPQB_H
#define WPQB_H

#include "common.h"

typedef union {
  alignas(EXTENDED_ALIGN_B) extended w;
  alignas(EXTENDED_ALIGN_B) struct {
    uint8_t a[10];
    /* p, q, b are stored in the highest, unused 2+2+2=6 bytes of w */
    uint16_t p, q, b;
  } i;
} wpqb;

extern wpqb *wpqb_alloc(const uint32_t n_a);
extern int wpqb_free(wpqb *const a);
extern int wpqb_init(wpqb *const a, const uint16_t p, const uint16_t q, const extended *const w);
extern int wpqb_print(FILE *const f, wpqb *const a);
extern int wpqb_dump(FILE *const f, wpqb *const a, const uint32_t n_a);
extern int wpqb_cmp(const wpqb *const a, const wpqb *const b);
extern int wpqb_sort(wpqb *const a, const uint32_t n_a);
extern int wpqb_ncp(wpqb *const a, const uint32_t n_a, uint16_t *const s, const uint16_t n_s, extended *const w);

#endif /* !WPQB_H */
