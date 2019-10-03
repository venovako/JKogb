#include "utils_.hpp"

#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_sort.h"

struct Qcmp {
  const qcmp _c;
  void *const _x;
  Qcmp(const qcmp c, void *const x) throw() : _c(c), _x(x) {}
  Qcmp(const Qcmp &q) throw() : _c(q._c), _x(q._x) {}
  bool operator()(const aw &a, const aw &b) const throw()
  {
    return
#ifdef _GNU_SOURCE
      (_c(&a, &b, _x) < 0)
#else // !_GNU_SOURCE
      (_c(_x, &a, &b) < 0)
#endif // ?_GNU_SOURCE
      ;
  }
};

int par_sort(const int t, aw *const a, const size_t n, const qcmp c, void *const x) throw()
{
  if (t < 0)
    return -1;
  if (!a)
    return -2;
  if (!n)
    return -3;
  if (!c)
    return -4;
  if (t)
    tbb::task_scheduler_init tsi(t);
  tbb::parallel_sort(a, (a + n), Qcmp(c, x));
  return 0;
}
