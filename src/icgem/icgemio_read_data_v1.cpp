#include "geodesy/core/geoconst.hpp"
#include "icgemio.hpp"
#include <cassert>
#include <charconv>

constexpr const int max_data_line = 256;

/** @warning coeffs should have already been initialized and allocated with
 *           enough memmory to hold the (to-be-) parsed coefficients.
 */
int dso::Icgem::parse_data_v1(int l, int k, const Icgem::Datetime &t,
                           dso::StokesCoeffs &coeffs) noexcept {

  /* check degree & order parameters */
  int error = 0;
  if (l > max_degree() || k > l) {
    fprintf(
        stderr,
        "[ERROR] Invalid degree/order given to data parse, i.e. (n,m)=(%d,%d) (traceback: %s)\n",
        l,k,__func__);
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
  double dt=0e0;
  int current_n=-1, current_m=-1;
  /* keep reading coefficients. Note that any entry of type 'trnd', 'asin' 
   * and 'acos' must be superceded by an entry of type 'gfct' so that we have 
   * the reference epoch 't0' (recorded in gfct) 
   */
  while (fin.getline(line, max_data_line) && !error) {
    error += resolve_icgem_data_line_v1(line, err, entry);
    {
      double C=0e0,S=0e0;
      if (entry.degree <= n && entry.order <= m) {
        /* set maximum degree/order */
        max_degree_collected = std::max(max_degree_collected, entry.degree);
        max_order_collected = std::max(max_order_collected, entry.order);
        /* compute C and S coeffs */
        switch (entry.key) {
          case DataEntryType::gfc:
            C = entry.C;
            S = entry.S;
            break;
          case DataEntryType::gfct:
            C = entry.C;
            S = entry.S;
            /* set current degree, order and t0 */
            current_n = entry.degree;
            current_m = entry.order;
            dt = t.diff<DateTimeDifferenceType::FractionalYears>(entry.t0);
            break;
          case DataEntryType::trnd:
            C = entry.C * dt;
            S = entry.S * dt;
            assert(current_n == entry.degree);
            assert(current_m == entry.order);
            break;
          case DataEntryType::asin:
            C = entry.C * std::sin((D2PI / entry.period) * dt);
            S = entry.S * std::sin((D2PI / entry.period) * dt);
            assert(current_n == entry.degree);
            assert(current_m == entry.order);
            break;
          case DataEntryType::acos:
            C = entry.C * std::cos((D2PI / entry.period) * dt);
            S = entry.S * std::cos((D2PI / entry.period) * dt);
            assert(current_n == entry.degree);
            assert(current_m == entry.order);
            break;
          default:
            assert(1==0);
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
