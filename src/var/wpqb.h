#ifndef WPQB_H
#define WPQB_H

#include <stdalign.h>
#include <stdint.h>

/* Intel 80-bit extended floating-point value stored in the lowest 10 bytes of a 16-byte variable */
typedef long double extended;

typedef union {
  alignas(16) extended w;
  alignas(16) struct {
    uint8_t a[10];
    /* p, q, b are stored in the highest, unused 2+2+2=6 bytes of w */
    uint16_t p, q, b;
  } i;
} wpqb;

extern int wpqb_init(wpqb *const a, const uint16_t p, const uint16_t q);
extern int wpqb_cmp(const wpqb *const a, const wpqb *const b);
extern int wpqb_sort(wpqb *const a, const uint32_t n_a);
extern int wpqb_ncp(wpqb *const a, const uint32_t n_a, uint16_t *const s, const uint16_t n_s);

#endif /* !WPQB_H */
