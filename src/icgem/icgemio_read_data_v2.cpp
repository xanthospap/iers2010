#include "geodesy/core/geoconst.hpp"
#include "icgemio.hpp"
#include <cassert>
#include <charconv>

constexpr const int max_data_line = 256;

/** @todo In this function, we recompute sin/cos values that are not needed 
 * (~ line 64 & 79) for C and S coefficients. This is branchless alright, but 
 * dould nevertheless. We could avoid this by checking the type of current 
 * coefficient parsed (i.e. 'trnd', 'asin', ...) and only compute trigs 
 * when needed (switch loop). See also the parse_data_v1 function.
 */
int dso::Icgem::parse_data_v2(int l, int k, const Icgem::Datetime &t,
                              dso::StokesCoeffs &coeffs) noexcept {

  /* check degree & order parameters */
  int error = 0;
  if (l > max_degree() || k > l) {
    fprintf(
        stderr,
        "[ERROR] Invalid degree/order given to data parse (traceback: %s)\n",
        __func__);
    error = 1;
  }

  /* set max degree and order we are collecting */
  const int n = l;
  const int m = k;
  
  /* clear out Stokes coeffs (i.e. set to zero) and resize (if needed) */
  coeffs.resize(n,m);
  coeffs.clear();

  /* max degree and order actually collected */
  int max_degree_collected = 0;
  int max_order_collected = 0;

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
      int interval_ok = (t >= entry.t0 && t < entry.t1);
      if (interval_ok && entry.degree <= n && entry.order <= m) {
        /* set maximum degree and order collected */
        max_degree_collected = std::max(max_degree_collected, entry.degree);
        max_order_collected = std::max(max_order_collected, entry.order);
        /* compute C and S coefficients */
        const double dt =
            t.diff<DateTimeDifferenceType::FractionalYears>(entry.t0);
        /* ------------------------------------------------------------------
         * WARNING !!
         * ------------------------------------------------------------------
         *  It took me fucking hours to debug this! if entry.period is zero, 
         *  the the computations below may result in a nan value, and cast the 
         *  collected values to garbage.
         *  Be extra careful to correctly avoid this issue.
         */
        double C = entry.C * (entry.key == DataEntryType::gfc) +
                   entry.C * (entry.key == DataEntryType::gfct) +
                   (entry.C * dt) * (entry.key == DataEntryType::trnd);
        if ((entry.key == DataEntryType::asin) ||
            (entry.key == DataEntryType::acos)) {
          C += (entry.C * std::sin((D2PI / entry.period) * dt)) *
                   (entry.key == DataEntryType::asin) +
               (entry.C * std::cos((D2PI / entry.period) * dt)) *
                   (entry.key == DataEntryType::acos);
        }
        double S = entry.S * (entry.key == DataEntryType::gfc) +
                   entry.S * (entry.key == DataEntryType::gfct) +
                   (entry.S * dt) * (entry.key == DataEntryType::trnd);
        if ((entry.key == DataEntryType::asin) ||
            (entry.key == DataEntryType::acos)) {
          S += (entry.S * std::sin((D2PI / entry.period) * dt)) *
                   (entry.key == DataEntryType::asin) +
               (entry.S * std::cos((D2PI / entry.period) * dt)) *
                   (entry.key == DataEntryType::acos);
        }
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

  /* no errors, assign further constants */
  coeffs.GM() = this->gm();
  coeffs.Re() = this->radius();
  coeffs.normalized() = this->is_normalized();
  /* if needed, set the actual dimensions of the StokesCoeffs instance */
  coeffs.shrink_dimensions(max_degree_collected, max_order_collected);

  /* all done */
  return 0;
}
