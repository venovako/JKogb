#include "tweight.h"

long double dtw(const double Aqp, const double Apq)
{
  const double Aqp_ = fabs(Aqp);
  const double Apq_ = fabs(Apq);
  double M, m;
  if (Aqp_ <= Apq_) {
    M = Apq_;
    m = Aqp_;
  }
  else {
    M = Aqp_;
    m = Apq_;
  }
  return ((M == 0.0) ? qNaN() : fmal(M, M, m * (long double)m));
}

long double ztw(const double _Complex Aqp, const double _Complex Apq)
{
  const double Aqp_r = fabs(creal(Aqp));
  const double Aqp_i = fabs(cimag(Aqp));
  const double Apq_r = fabs(creal(Apq));
  const double Apq_i = fabs(cimag(Apq));
  double M_r, m_r, M_i, m_i, t;
  if (Aqp_r <= Apq_r) {
    M_r = Apq_r;
    m_r = Aqp_r;
  }
  else {
    M_r = Aqp_r;
    m_r = Apq_r;
  }
  if (Aqp_i <= Apq_i) {
    M_i = Apq_i;
    m_i = Aqp_i;
  }
  else {
    M_i = Aqp_i;
    m_i = Apq_i;
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
  // now, M_r >= M_i >= m_r >= m_i
  return ((M_r == 0.0) ? qNaN() : fmal(M_r, M_r, fmal(M_i, M_i, fmal(m_r, m_r, m_i * (long double)m_i))));
}
