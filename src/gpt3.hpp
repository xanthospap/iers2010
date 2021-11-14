#ifndef __GPT3_TROPO_HPP__
#define __GPT3_TROPO_HPP__

#include "ggdatetime/dtcalendar.hpp"

namespace dso {
namespace gpt3 {
constexpr int GPT3_5_GRID_LINES = 2593; // including header line

struct gpt3_5_grid {
  double p_grid[GPT3_5_GRID_LINES][5];
  double T_grid[GPT3_5_GRID_LINES][5];
  double Q_grid[GPT3_5_GRID_LINES][5];
  double dT_grid[GPT3_5_GRID_LINES][5];
  double u_grid[GPT3_5_GRID_LINES];
  double Hs_grid[GPT3_5_GRID_LINES];
  double ah_grid[GPT3_5_GRID_LINES][5];
  double aw_grid[GPT3_5_GRID_LINES][5];
  double la_grid[GPT3_5_GRID_LINES][5];
  double Tm_grid[GPT3_5_GRID_LINES][5];
  double Gn_h_grid[GPT3_5_GRID_LINES][5];
  double Ge_h_grid[GPT3_5_GRID_LINES][5];
  double Gn_w_grid[GPT3_5_GRID_LINES][5];
  double Ge_w_grid[GPT3_5_GRID_LINES][5];
}; // gpt3_5_grid

int parse_gpt3_5_grid(const char *gridfn, gpt3_5_grid *grid) noexcept;
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
};

int gpt3_5_fast(const dso::datetime<dso::nanoseconds> &t, const double *lat,
                const double *lon, const double *hell, int num_stations, int it,
                const char *grid_file, gpt3_result *g3out) noexcept;

} // namespace dso

#endif
