#include "wpqb.h"

static_assert(sizeof(wpqb) == sizeof(long double), "sizeof(wpqb) != sizeof(long double)");

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
  if (a->w == a->w)
    return -5;
  // both weights are NaN
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

#ifdef _OPENMP
static void wpqb_sort1(const uint32_t n_a, wpqb a[static 1])
{
  if (!n_a)
    return;
  const uint32_t na = (n_a - 1u);

  for (uint32_t swps = ~0u, s = 0u, oe = 0u; swps; (oe ^= 1u), (swps = s), (s = 0u)) {
    // odd-even sort
#pragma omp parallel for default(none) shared(oe,na,a) reduction(+:s)
    for (uint32_t i = oe; i < na; i += 2u) {
      const uint32_t j = (i + 1u);
      if (wpqb_cmp((a + i), (a + j)) > 0) {
        const wpqb t = a[i];
        a[i] = a[j];
        a[j] = t;
        ++s;
      }
    }
  }
}
#else /* !_OPENMP */
static void wpqb_sort0(const uint32_t n_a, wpqb a[static 1])
{
  if (n_a)
    qsort(a, n_a, sizeof(wpqb), (int (*)(const void*, const void*))wpqb_cmp);
}
#endif /* ?_OPENMP */

static void wpqb_ncp0(const uint16_t n, const uint32_t n_a, wpqb a[static 1], const uint16_t n_s, uint32_t s[static 1], wpqb_info w[static 1])
{
  w->w = qNaN();
  w->i.s = UINT16_C(0);
  w->i.f = 0u;

  if (!n)
    return;
  if (!n_s)
    return;
  if (!n_a)
    return;

  bool r[n];
  for (uint16_t i = UINT16_C(0); i < n; ++i)
    r[i] = false;

  uint16_t ns = (uint16_t)(n >> 1u);
  ns = (uint16_t)(((uint32_t)ns <= n_a) ? ns : (uint16_t)n_a);
  ns = (uint16_t)((ns <= n_s) ? ns : n_s);

  for (uint32_t j = 0u, k; w->i.s < ns; j = k) {
    s[w->i.s] = j;
    r[a[j].i.p] = true;
    r[a[j].i.q] = true;
    for (k = j + 1u; k < n_a; ++k)
      if (!(r[a[k].i.p] || r[a[k].i.q]))
        break;
    (w->i.s)++;
    if (k >= n_a)
      break;
  }

  w->w = 0.0L;
  for (uint16_t i = w->i.s; i; )
    w->w += a[s[--i]].w;
}

#ifdef _OPENMP
static void wpqb_update(uint32_t s[static restrict 1], const uint32_t my_s[static restrict 1], wpqb_info w[static restrict 1], const wpqb_info my_w[static restrict 1])
{
  if (!(my_w->w <= w->w)) {
    for (uint16_t i = UINT16_C(0); i < my_w->i.s; ++i)
      s[i] = my_s[i] + my_w->i.f;
    *w = *my_w;
  }
}

static void wpqb_ncpt(const uint16_t n, const uint32_t n_a, wpqb a[static 1], const uint32_t f, const uint16_t n_s, uint32_t s[static 1], wpqb_info w[static 1])
{
  uint32_t my_s[n_s];
  wpqb_info my_w;
  wpqb_ncp0(n, n_a, a, n_s, my_s, &my_w);
  if (my_w.i.s) {
    my_w.i.f = f;
#pragma omp critical
    wpqb_update(s, my_s, w, &my_w);
  }
}

static void wpqb_ncp1(const uint16_t n, const uint32_t n_a, wpqb a[static 1], const uint16_t n_s, uint32_t s[static 1], wpqb_info w[static 1])
{
  w->w = qNaN();
  w->i.s = UINT16_C(0);
  w->i.f = 0u;

  if (!n)
    return;
  if (!n_s)
    return;
  if (!n_a)
    return;

#ifdef __ICL
#pragma omp parallel for default(none) shared(n,n_a,a,n_s,s,w) schedule(dynamic,1)
#else /* !__ICL */
#pragma omp parallel for default(none) shared(n,n_a,a,n_s,s,w) schedule(nonmonotonic:dynamic,1)
#endif /* ?__ICL */
  for (uint32_t i = 0u; i < n_a; ++i)
    wpqb_ncpt(n, n_a - i, a + i, i, n_s, s, w);
}
#endif /* _OPENMP */

void wpqb_ncp(const uint16_t n, const uint32_t n_a, wpqb a[static 1], const uint16_t n_s, uint32_t s[static 1], wpqb_info w[static 1])
{
  w->w = qNaN();
  w->i.s = UINT16_C(0);
  w->i.f = 0u;

#ifdef _OPENMP
  wpqb_sort1(n_a, a);
#else /* !_OPENMP */
  wpqb_sort0(n_a, a);
#endif /* ?_OPENMP */

  for (w->i.f = n_a; w->i.f; ) {
    --(w->i.f);
    if (a[w->i.f].w >= -LDBL_MAX) {
      ++(w->i.f);
      break;
    }
  }
  if (!(w->i.f))
    return;

#ifdef _OPENMP
  wpqb_ncp1(n, w->i.f, a, n_s, s, w);
#else /* !_OPENMP */
  wpqb_ncp0(n, w->i.f, a, n_s, s, w);
#endif /* ?_OPENMP */
}
