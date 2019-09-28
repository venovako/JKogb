#include "utils_.hpp"

#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_sort.h"

struct Qcmp {
  const qcmp _c;
  Qcmp(const qcmp c) throw() : _c(c) {}
  Qcmp(const Qcmp &q) throw() : _c(q._c) {}
  bool operator()(const aw &a, const aw &b) const throw() { return (_c(&a, &b) < 0); }
};

void par_sort(const int t, aw *const a, const size_t n, const qcmp c) throw()
{
  if (!a)
    return;
  if (!n)
    return;
  if (!c)
    return;
  tbb::task_scheduler_init tsi(t);
  tbb::parallel_sort(a, a+n, Qcmp(c));
}
