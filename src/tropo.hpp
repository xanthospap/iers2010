#ifndef __TROPO_MODELS_ALTERNATIVES_HPP__
#define __TROPO_MODELS_ALTERNATIVES_HPP__

#include "ggdatetime/dtcalendar.hpp"
#ifdef DEBUG
#include <cstdio>
#endif

namespace dso {
double asknewet(double e, double Tm, double lambda) noexcept;
double saasthyd(double p, double dlat, double hell) noexcept;
int vmf3(double ah, double aw, dso::datetime<dso::nanoseconds> &t, double lat,
         double lon, double zd, double &mfh, double &mfw) noexcept;
namespace gpt3 {

enum class Gpt3Grid : char { grid1x1, grid5x5 };
template <Gpt3Grid G> struct gpt3_grid_attributes {
};

template <> struct gpt3_grid_attributes<Gpt3Grid::grid1x1> {
  template <typename T> static constexpr T grid_tick() { return T(1e0); }
  static constexpr unsigned num_lines = 64801-1;
};
template <> struct gpt3_grid_attributes<Gpt3Grid::grid5x5> {
  template <typename T> static constexpr T grid_tick() { return T(5e0); }
  static constexpr unsigned num_lines = 2593-1;
};


struct gpt3_grid {

  gpt3_grid(int rows=0) {
    if (rows) allocate(rows);
  };
  ~gpt3_grid() noexcept {
    if (size) dealloc();
  }
  gpt3_grid(const gpt3_grid&) noexcept = delete;
  gpt3_grid(gpt3_grid&&) noexcept = delete;
  gpt3_grid& operator=(const gpt3_grid&) noexcept = delete;
  gpt3_grid& operator=(gpt3_grid&&) noexcept = delete;

  void allocate(unsigned num_rows);

  void dealloc() noexcept;

  unsigned size = 0;
  double **p_grid=nullptr;
  double **T_grid=nullptr;
  double **Q_grid=nullptr;
  double **dT_grid=nullptr;
  double *u_grid=nullptr;
  double *Hs_grid=nullptr;
  double **ah_grid=nullptr;
  double **aw_grid=nullptr;
  double **la_grid=nullptr;
  double **Tm_grid=nullptr;
  double **Gn_h_grid=nullptr;
  double **Ge_h_grid=nullptr;
  double **Gn_w_grid=nullptr;
  double **Ge_w_grid=nullptr;
};

int parse_gpt3_grid(const char *gridfn, gpt3_grid *grid) noexcept;
} // namespace gpt3

struct gpt3_result {
  double p;    // pressure in hPa
  double T;    // temperature in degrees Celsius
  double dT;   // temperature lapse rate in degrees per km
  double Tm;   // mean temperature weighted with the water vapor in degrees
               // Kelvin
  double e;    // water vapour pressure in hPa
  double ah;   // hydrostatic mapping function coefficient(VMF3)
  double aw;   // wet mapping function coefficient(VMF3)
  double la;   // water vapour decrease factor
  double undu; // geoid undulation in m
  double Gn_h; // hydrostatic north gradient in m
  double Ge_h; // hydrostatic east gradient in m
  double Gn_w; // wet north gradient in m
  double Ge_w; // wet east gradient in m
};             // gpt3_result

struct vmf3_hw {
  double mfh, mfw;
}; // vmf3_result

int gpt3_fast(const dso::datetime<dso::nanoseconds> &t, const double *lat,
                const double *lon, const double *hell, int num_stations, int it,
                const gpt3::gpt3_grid *grid5x5, gpt3_result *g3out) noexcept;
int gpt3_fast(const dso::datetime<dso::nanoseconds> &t, const double *lat,
                const double *lon, const double *hell, int num_stations, int it,
                const char *grid_file, gpt3_result *g3out, int &grid_step) noexcept;
int vmf3(const dso::gpt3_result *gptres, dso::datetime<dso::nanoseconds> &t,
         const double *vlat, const double *vlon, const double *vzd,
         dso::vmf3_hw *vmfres, int num_sta) noexcept;

} // namespace dso
#endif
