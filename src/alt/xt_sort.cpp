#include "xt_sort.h"
#include <algorithm>
#include <oneapi/tbb/parallel_sort.h>

struct xt {
  char x[10];
  char p[3];
  char q[3];
};

static_assert(sizeof(xt) == sizeof(long double));

static bool cmp(const long double &a, const long double &b) throw()
{
  const long double *const pa = &a;
  const long double *const pb = &b;
  if (pa == pb)
    return false;
  if (a > b)
    return true;
  if (a == b) {
    const xt &xa = *(const xt*)pa;
    const xt &xb = *(const xt*)pb;
    unsigned ap = 0u, aq = 0u, bp = 0u, bq = 0u;
    char *c = (char*)&ap;
    c[0] = xa.p[0];
    c[1] = xa.p[1];
    c[2] = xa.p[2];
    c = (char*)&aq;
    c[0] = xa.q[0];
    c[1] = xa.q[1];
    c[2] = xa.q[2];
    c = (char*)&bp;
    c[0] = xb.p[0];
    c[1] = xb.p[1];
    c[2] = xb.p[2];
    c = (char*)&bq;
    c[0] = xb.q[0];
    c[1] = xb.q[1];
    c[2] = xb.q[2];
    const unsigned ab = aq - ap;
    const unsigned bb = bq - bp;
    if (ab > bb)
      return true;
    if (ab == bb)
      return (aq > bq);
  }
  return false;
}

void xt_sort_(const MKL_INT *const n, long double *const x) throw()
{
  if (x) {
    if (*n > 0)
      oneapi::tbb::parallel_sort(x, (x + *n), cmp);
    else if (*n < 0)
      std::sort(x, (x - *n), cmp);
  }
}
