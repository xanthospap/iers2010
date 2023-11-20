#include "geodesy/geoconst.hpp"
#include "icgemio.hpp"
#include <cassert>
#include <charconv>

constexpr const int max_data_line = 256;

/** @warning coeffs should have already been initialized and allocated with
 *           enough memmory to hold the (to-be-) parsed coefficients.
 */
int dso::Icgem::parse_data(int l, int k, const Icgem::Datetime &t,
                           dso::StokesCoeffs &coeffs) noexcept {

  /* clear out Stokes coeffs (i.e. set to zero) */
  coeffs.clear();

  /* check degree & order parameters */
  int error = 0;
  if (l > max_degree() || k > l) {
    fprintf(
        stderr,
        "[ERROR] Invalid degree/order given to data parse (traceback: %s)\n",
        __func__);
    error = 1;
  }
  if (coeffs.max_degree() < l) {
    fprintf(stderr,
            "[ERROR] Cannot read harmonics of degree %d to HarmonicsCoeffs of "
            "degree %d (traceback: %s)\n",
            l, coeffs.max_degree(), __func__);
    error = 1;
  }
  if (coeffs.max_order() < k) {
    fprintf(stderr,
            "[ERROR] Cannot read harmonics of order %d to HarmonicsCoeffs of "
            "order %d (traceback: %s)\n",
            k, coeffs.max_order(), __func__);
    error = 1;
  }

  /* set max degree and order we are collecting */
  const int n = l;
  const int m = k;

  /* go the start of the data section (in file) -- assert file is open */
  std::ifstream fin(filename().c_str());
  if (!fin.is_open()) {
    fprintf(stderr, "[ERROR] Failed opening icgem file %s (traceback: %s)\n",
            filename().c_str(), __func__);
    error = 1;
  }

  /* stop in case of error */
  if (error)
    return error;

  fin.seekg(data_section_pos);
  char line[max_data_line];
  Icgem::DataEntry entry;
  const auto err = this->errors();
  /* keep reading coefficients */
  while (fin.getline(line, max_data_line) && !error) {
    error += resolve_icgem_data_line_v2(line, err, entry);
    {
      int interval_ok = (t > entry.t0 && t <= entry.t1);
      if (interval_ok && entry.degree <= n && entry.order <= m) {
        const double dt =
            t.diff<DateTimeDifferenceType::FractionalYears>(entry.t0);
        const double C = entry.C * (entry.key == DataEntryType::gfc) +
                         (entry.C * dt) * (entry.key == DataEntryType::gfct) +
                         (entry.C * std::sin((D2PI / entry.period) * dt)) *
                             (entry.key == DataEntryType::asin) +
                         (entry.C * std::cos((D2PI / entry.period) * dt)) *
                             (entry.key == DataEntryType::acos);
        const double S = entry.C * (entry.key == DataEntryType::gfc) +
                         (entry.S * dt) * (entry.key == DataEntryType::gfct) +
                         (entry.S * std::sin((D2PI / entry.period) * dt)) *
                             (entry.key == DataEntryType::asin) +
                         (entry.S * std::cos((D2PI / entry.period) * dt)) *
                             (entry.key == DataEntryType::acos);
        coeffs.C(entry.degree, entry.order) += C;
        coeffs.S(entry.degree, entry.order) += S;
      } else {
        ;
      }
    }
  }

  /* check for errors */
  if (error) {
    fprintf(
        stderr,
        "[ERROR] Failed parsing coefficients from file %s (traceback: %s)\n",
        filename().c_str(), __func__);
    fprintf(stderr, "[ERROR] Last line read was [%s] (traceback: %s)\n", line,
            __func__);
    return 1;
  }
  if (!fin) {
    if (fin.eof())
      fin.clear();
    else {
      fprintf(stderr,
              "[ERROR] Failed parsing coefficients from file %s (traceback: "
              "%s)\n",
              filename().c_str(), __func__);
      fprintf(stderr, "[ERROR] Some stream error occured (traceback: %s)\n",
              __func__);
      return 1;
    }
  }

  /* all done */
  return 0;
}
