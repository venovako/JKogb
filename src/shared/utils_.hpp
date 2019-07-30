#ifndef UTILS_HPP
#define UTILS_HPP

#include <cstdint>
#include <cstdlib>

struct aw {
  double w;
  int64_t p, q, b;
};

typedef int (*qcmp)(const aw *const, const aw *const);
extern "C" void par_sort(aw *const a, const size_t n, const qcmp c) throw();

#endif // !UTILS_HPP
