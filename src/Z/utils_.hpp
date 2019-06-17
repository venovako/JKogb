#ifndef UTILS_HPP
#define UTILS_HPP

#include <cstdint>
#include <cstdlib>

struct dzbw {
  double w;
  int64_t p, q, b;
};

typedef int (*qcmp)(const dzbw *const, const dzbw *const);
extern "C" void par_sort(dzbw *const a, const size_t n, const qcmp c);

#endif // !UTILS_HPP
