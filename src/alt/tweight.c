#include "tweight.h"

long double dfn2(const double a, const double b)
{
  const double a_ = fabs(a);
  const double b_ = fabs(b);
  double M, m;
  if (a_ <= b_) {
    M = b_;
    m = a_;
  }
  else {
    M = a_;
    m = b_;
  }
  if (m == 0.0) {
    if (M == 0.0)
      return 0.0L;
    return (m * (long double)m);
  }
  return fmal(M, M, m * (long double)m);
}

long double zfn2(const double _Complex a, const double _Complex b)
{
  const double a_i = fabs(cimag(a));
  const double b_i = fabs(cimag(b));
  const double a_r = fabs(creal(a));
  const double b_r = fabs(creal(b));
  double M_r, m_r, M_i, m_i, t;
  if (a_r <= b_r) {
    M_r = b_r;
    m_r = a_r;
  }
  else {
    M_r = a_r;
    m_r = b_r;
  }
  if (a_i <= b_i) {
    M_i = b_i;
    m_i = a_i;
  }
  else {
    M_i = a_i;
    m_i = b_i;
  }
  if (M_r < M_i) {
    t = M_r;
    M_r = M_i;
    M_i = t;
  }
  if (m_r < m_i) {
    t = m_r;
    m_r = m_i;
    m_i = t;
  }
  if (M_i < m_r) {
    t = M_i;
    M_i = m_r;
    m_r = t;
    if (m_r < m_i) {
      t = m_r;
      m_r = m_i;
      m_i = t;
    }
  }
  // now, M_r >= M_i >= m_r >= m_i >= 0
  if (m_i == 0.0) {
    if (m_r == 0.0) {
      if (M_i == 0.0) {
        if (M_r == 0.0)
          return 0.0L;
        return (M_r * (long double)M_r);
      }
      return fmal(M_r, M_r, M_i * (long double)M_i);
    }
    return fmal(M_r, M_r, fmal(M_i, M_i, m_r * (long double)m_r));
  }
  return fmal(M_r, M_r, fmal(M_i, M_i, fmal(m_r, m_r, m_i * (long double)m_i)));
}

long double dtw(const double Aqp, const double Apq)
{
  const long double w = dfn2(Aqp, Apq);
  return ((w == 0.0L) ? qNaN() : w);
}

long double ztw(const double _Complex Aqp, const double _Complex Apq)
{
  const long double w = zfn2(Aqp, Apq);
  return ((w == 0.0L) ? qNaN() : w);
}
