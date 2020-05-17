#include "wpqb.h"

int wpqb_cmp(const wpqb a[static 1], const wpqb b[static 1])
{
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

uint32_t wpqb_clean(const uint32_t n_a, wpqb a[static 1])
{
  if (!n_a)
    return 0u;

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

void wpqb_sort0(const uint32_t n_a, wpqb a[static 1])
{
  qsort(a, n_a, sizeof(wpqb), (int (*)(const void*, const void*))wpqb_cmp);
}

void wpqb_ncp0(const uint16_t n, const uint32_t n_a, wpqb a[static 1], const uint16_t n_s, uint32_t s[static 1], wpqb_info w[static 1])
{
  for (unsigned i = 0u; i < 10u; ++i)
    (w->i.a)[i] = UINT8_C(0xFF);
  w->i.s = UINT16_C(0);
  w->i.f = 0u;

  if (!n)
    return;
  if (!n_s)
    return;
  if (!n_a)
    return;

  bool p[n];
  for (uint16_t i = UINT16_C(0); i < n; ++i)
    p[i] = false;

  bool q[n];
  for (uint16_t i = UINT16_C(0); i < n; ++i)
    q[i] = false;

  uint16_t ns = (n >> 1u);
  ns = ((ns <= n_a) ? ns : (uint16_t)n_a);
  ns = ((ns <= n_s) ? ns : n_s);

  for (uint32_t j = 0u, k; w->i.s < ns; j = k) {
    s[w->i.s++] = j;
    q[a[j].i.q] = p[a[j].i.p] = true;
    for (k = j + 1u; k < n_a; ++k)
      if (!(p[a[k].i.p] || q[a[k].i.q]))
        break;
    if (k >= n_a)
      break;
  }

  w->w = 0.0L;
  for (uint16_t i = w->i.s; i; )
    w->w += a[s[--i]].w;
}

static void wpqb_update(uint32_t s[static 1], const uint32_t my_s[static 1], wpqb_info w[static 1], const wpqb_info my_w[static 1])
{
  if (!(my_w->w <= w->w)) {
    for (uint16_t i = UINT16_C(0); i < my_w->i.s; ++i)
      s[i] = my_s[i];
    *w = *my_w;
  }
}

static void wpqb_ncpt(const uint16_t n, const uint32_t n_a, wpqb a[static 1], const uint32_t f, const uint16_t n_s, uint32_t s[static 1], wpqb_info w[static 1])
{
  wpqb_info my_w;
  for (unsigned i = 0u; i < 10u; ++i)
    my_w.i.a[i] = UINT8_C(0xFF);
  my_w.i.s = UINT16_C(0);
  my_w.i.f = f;

  uint32_t my_s[n_s];
  wpqb_ncp0(n, n_a, a, n_s, my_s, &my_w);
#pragma omp critical
  wpqb_update(s, my_s, w, &my_w);
}

void wpqb_ncp1(const uint16_t n, const uint32_t n_a, wpqb a[static 1], const uint16_t n_s, uint32_t s[static 1], wpqb_info w[static 1])
{
  for (unsigned i = 0u; i < 10u; ++i)
    (w->i.a)[i] = UINT8_C(0xFF);
  w->i.s = UINT16_C(0);
  w->i.f = 0u;

  if (!n)
    return;
  if (!n_s)
    return;
  if (!n_a)
    return;

#ifdef OLD_OMP
#pragma omp parallel for default(none) shared(n,n_a,a,n_s,s,w) schedule(dynamic,1)
#else /* !OLD_OMP */
#pragma omp parallel for default(none) shared(n,n_a,a,n_s,s,w) schedule(nonmonotonic:dynamic,1)
#endif /* ?OLD_OMP */
  for (uint32_t i = 0u; i < n_a; ++i)
    wpqb_ncpt(n, n_a - i, a + i, i, n_s, s, w);
}
