#include "stokes_coefficients.hpp"
#include <stdexcept>

void dso::StokesCoeffs::resize(int degree, int order) noexcept {
  _Cnm.resize(degree + 1, degree + 1);
  _Snm.resize(degree + 1, degree + 1);
  m_degree = degree;
  m_order = order;
}

void dso::StokesCoeffs::cresize(int degree, int order) noexcept {
  /* first copy the data to a new instance */
  _Cnm.cresize(degree + 1, degree + 1);
  _Snm.cresize(degree + 1, degree + 1);
  m_degree = degree;
  m_order = order;
}

dso::StokesCoeffs &dso::StokesCoeffs::operator+=(const dso::StokesCoeffs &sc) {
  if (m_degree < sc.max_degree() || m_order < sc.max_order()) {
    throw std::runtime_error("[ERROR] Failed to add Stokes coefficients\n");
  }

  for (int m = 0; m <= sc.max_order(); m++) {
    const double *colCnm2 = sc.Cnm().column(m);
    double *colCnm = this->Cnm().column(m);
    for (int n = m; n <= sc.max_degree(); n++) {
      colCnm[n] += colCnm2[n];
    }
  }

  for (int m = 0; m <= sc.max_order(); m++) {
    const double *colSnm2 = sc.Snm().column(m);
    double *colSnm = this->Snm().column(m);
    for (int n = m; n <= sc.max_degree(); n++) {
      colSnm[n] += colSnm2[n];
    }
  }

  return *this;
}
