#include "utils_.hpp"

#include "tbb/parallel_sort.h"

struct Qcmp {
  const qcmp _c;
  Qcmp(const qcmp c) throw() : _c(c) {}
  Qcmp(const Qcmp &q) throw() : _c(q._c) {}
  bool operator()(const aw &a, const aw &b) const throw() { return (_c(&a, &b) < 0); }
};

void par_sort(aw *const a, const size_t n, const qcmp c) throw()
{
  if (!a)
    return;
  if (!n)
    return;
  if (!c)
    return;
  tbb::parallel_sort(a, a+n, Qcmp(c));
}
