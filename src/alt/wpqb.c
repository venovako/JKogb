#include "wpqb.h"

int wpqb_cmp(const wpqb a[static 1], const wpqb b[static 1])
{
  assert(a);
  assert(b);
  if (a == b)
    return 0;
  if (a->w < b->w)
    return 1;
  if (a->w > b->w)
    return -1;
  if (a->w == b->w) {
    if (a->i.b < b->i.b)
      return 2;
    if (a->i.b > b->i.b)
      return -2;
    if (a->i.p < b->i.p)
      return 3;
    if (a->i.p > b->i.p)
      return -3;
    if (a->i.q < b->i.q)
      return 4;
    if (a->i.q > b->i.q)
      return -4;
    return 0;
  }
  if (b->w == b->w)
    return 5;
  return -5;
}

uint32_t wpqb_clean(wpqb *const a, const uint32_t n_a)
{
  if (!n_a)
    return 0u;
  if (!a)
    return (n_a + 1u);
  // ii: index of the last non-invalid
  uint32_t ii = n_a - 1u;
  while (wpqb_invalid(a + ii)) {
    if (ii)
      --ii;
    else
      return 0u;
  }
  uint32_t nn = ii + 1u;
  for (uint32_t i = 0u; i < ii; ++i) {
    if (wpqb_invalid(a + i)) {
      const wpqb t = a[i];
      a[i] = a[ii];
      a[ii] = t;
      --nn;
      for (--ii; !(a[ii].w == a[ii].w); --ii) /**/;
    }
  }
  return nn;
}

void wpqb_ncp0(const uint16_t n, const uint32_t n_a, wpqb *const restrict a, const uint32_t f, const uint16_t n_s, uint32_t *const restrict s, wpqb_info *const restrict w)
{
  assert(w);
  for (unsigned i = 0u; i < 10u; ++i)
    (w->i.a)[i] = UINT8_C(0xFF);
  w->i.s = UINT16_C(0);
  w->i.f = f;

  if (!n)
    return;
  if (!n_s)
    return;
  if (!n_a)
    return;
  if (f >= n_a)
    return;
  assert(a);
  assert(s);

  bool p[n];
  for (uint16_t i = UINT16_C(0); i < n; ++i)
    p[i] = false;

  bool q[n];
  for (uint16_t i = UINT16_C(0); i < n; ++i)
    q[i] = false;

  uint16_t ns = (n >> 1u);
  ns = ((ns <= n_a) ? ns : (uint16_t)n_a);
  ns = ((ns <= n_s) ? ns : n_s);

  w->w = 0.0L;
  for (uint32_t j = f, k; w->i.s < ns; j = k) {
    w->w += a[s[(w->i.s)++] = j].w;
    q[a[j].i.q] = p[a[j].i.p] = true;
    for (k = j + 1u; k < n_a; ++k)
      if (!(p[a[k].i.p] || q[a[k].i.q]))
        break;
    if (k >= n_a)
      break;
  }
}

static void wpqb_update(uint32_t *const restrict s, const uint32_t *const restrict my_s, wpqb_info *const restrict w, const wpqb_info *const restrict my_w)
{
  if (!(my_w->w <= w->w)) {
    for (uint16_t i = UINT16_C(0); i < my_w->i.s; ++i)
      s[i] = my_s[i];
    *w = *my_w;
  }
}

static void wpqb_ncp(const uint16_t n, const uint32_t n_a, wpqb *const restrict a, const uint32_t f, const uint16_t n_s, uint32_t *const restrict s, wpqb_info *const restrict w)
{
  wpqb_info my_w;
  uint32_t *const my_s = (uint32_t*)alloca(n_s * sizeof(uint32_t));
  assert(my_s);
  wpqb_ncp0(n, n_a, a, f, n_s, my_s, &my_w);
#pragma omp critical
  wpqb_update(s, my_s, w, &my_w);
}

void wpqb_ncp1(const uint16_t n, const uint32_t n_a, wpqb *const restrict a, const uint32_t f, const uint16_t n_s, uint32_t *const restrict s, wpqb_info *const restrict w)
{
  assert(w);
  for (unsigned i = 0u; i < 10u; ++i)
    (w->i.a)[i] = UINT8_C(0xFF);
  w->i.s = UINT16_C(0);
  w->i.f = f;

  if (!n)
    return;
  if (!n_s)
    return;
  if (!n_a)
    return;
  if (f >= n_a)
    return;
  assert(a);
  assert(s);

  const uint32_t na = n_a - f;
#pragma omp parallel for default(none) shared(f,na,n,n_a,a,n_s,s,w) schedule(dynamic,1)
  for (uint32_t i = f; i < na; ++i)
    wpqb_ncp(n, n_a, a, i, n_s, s, w);
}
