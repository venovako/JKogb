#ifndef UTILS_HPP
#define UTILS_HPP

#include <cstdint>
#include <cstdlib>

struct aw {
  double w;
  int64_t p, q, b;
};

#ifdef _GNU_SOURCE
typedef int (*qcmp)(const aw *const, const aw *const, void *const);
#else // !_GNU_SOURCE
typedef int (*qcmp)(void *const, const aw *const, const aw *const);
#endif // ?_GNU_SOURCE
extern "C" int par_sort(const int t, aw *const a, const size_t n, const qcmp c, void *const x) throw();

#endif // !UTILS_HPP
