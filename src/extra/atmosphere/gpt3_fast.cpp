#include "datetime/dtcalendar.hpp"
#include "geodesy/units.hpp"
#include "iersc.hpp"
#include "tropo.hpp"
#include <cmath>
#include <fstream>
#ifdef DEBUG
#include <cassert>
#endif

using dso::gpt3::gpt3_grid_attributes;
using dso::gpt3::Gpt3GridResolution;

template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

int dso::gpt3_fast(double fractional_doy,
                   const std::vector<std::array<double, 3>> &ellipsoidal,
                   int it, const dso::Gpt3Grid &grid,
                   std::vector<dso::gpt3_result> &g3out) noexcept {

  g3out.clear();
  g3out = std::vector<dso::gpt3_result>(ellipsoidal.size(), dso::gpt3_result());

  // grid info
  double grid_tick, half_grid_tick;
  int ilon_max, ipod_max;
  if (grid.num_rows ==
      gpt3_grid_attributes<Gpt3GridResolution::grid1x1>::num_lines) {
    grid_tick = gpt3_grid_attributes<
        Gpt3GridResolution::grid1x1>::template grid_tick<double>();
    half_grid_tick =
        gpt3_grid_attributes<Gpt3GridResolution::grid1x1>::template grid_tick<
            double>() /
        2e0;
    ipod_max = 180;
    ilon_max = 360;
  } else if (grid.num_rows ==
             gpt3_grid_attributes<Gpt3GridResolution::grid5x5>::num_lines) {
    grid_tick = gpt3_grid_attributes<
        Gpt3GridResolution::grid5x5>::template grid_tick<double>();
    half_grid_tick =
        gpt3_grid_attributes<Gpt3GridResolution::grid5x5>::template grid_tick<
            double>() /
        2e0;
    ipod_max = 36;
    ilon_max = 72;
  } else {
    fprintf(stderr,
            "[ERROR] gpt3 grid file seems to have invalid number of rows! "
            "(traceback: %s)\n",
            __func__);
    return 15;
  }

  // determine the GPT3 coefficients
  // ------------------------------------------------------------------------
  constexpr double pi(iers2010::DPI);
  // mean gravity in m/s**2
  constexpr double gm = 9.80665e0;
  // molar mass of dry air in kg/mol
  constexpr double dMtr = 28.965e-3;
  // universal gas constant in J/K/mol
  constexpr double Rg = 8.3143e0;
  const double fdoy = fractional_doy;

  double cosfy = 0e0, coshy = 0e0, sinfy = 0e0, sinhy = 0e0;
  // factors for amplitudes
  if (it != 1) {
    cosfy = std::cos(fdoy / 365.25e0 * 2 * pi); // coefficient for A1
    coshy = std::cos(fdoy / 365.25e0 * 4 * pi); // coefficient for B1
    sinfy = std::sin(fdoy / 365.25e0 * 2 * pi); // coefficient for A2
    sinhy = std::sin(fdoy / 365.25e0 * 4 * pi); // coefficient for B2
  }

  // used for indexing lines
  int indx[4];

  // loop over stations
  for (auto ell = ellipsoidal.cbegin(); ell != ellipsoidal.cend(); ++ell) {
    int k = std::distance(ellipsoidal.cbegin(), ell);
    const double lon = ell->operator[](0);
    const double lat = ell->operator[](1);
    const double hell = ell->operator[](2);

    // only positive longitude in degrees
    const double plon = dso::rad2deg(lon + (lon < 0) * 2e0 * pi);
    // transform to polar distance in degrees
    const double ppod = dso::rad2deg(-lat + pi / 2e0);

    // find the index (line in the grid file) of the nearest point
    int ipod = std::floor((ppod + grid_tick) / grid_tick);
    int ilon = std::floor((plon + grid_tick) / grid_tick);

    // normalized (to one) differences, can be positive or negative
    const double diffpod =
        (ppod - (ipod * grid_tick - half_grid_tick)) / grid_tick;
    const double difflon =
        (plon - (ilon * grid_tick - half_grid_tick)) / grid_tick;
    if (ipod > ipod_max)
      ipod = ipod_max;
    if (ilon > ilon_max)
      ilon = 1;
    else if (ilon == 0)
      ilon = ilon_max;

    // get the number of the corresponding line
    indx[0] = (ipod - 1) * ilon_max + ilon - 1;

    // near the poles: nearest neighbour interpolation, otherwise: bilinear
    int bilinear =
        (ppod > half_grid_tick && ppod < 180e0 - half_grid_tick) ? 1 : 0;

    // case of nearest neighbourhood
    if (!bilinear) {
      int ix = indx[0];
#ifdef DEBUG
      if (ix < 0 || ix >= (int)grid.num_rows) {
        fprintf(stderr,
                ">> error! Index=%d and size is %u; ipod=%d, ilon=%d, "
                "plon=%.5f, ppod=%.5f\n",
                ix, grid.num_rows, ipod, ilon, plon, ppod);
      }
      assert(ix >= 0 && ix < (int)grid.num_rows);
#endif

      // transforming ellipsoidal height to orthometric height
      g3out[k].undu = grid.u_grid(ix);
      const double hgt = hell - g3out[k].undu;

      // pressure, temperature at the height of the grid
      const double T0 = grid.t_grid(ix, 0) + grid.t_grid(ix, 1) * cosfy +
                        grid.t_grid(ix, 2) * sinfy +
                        grid.t_grid(ix, 3) * coshy + grid.t_grid(ix, 4) * sinhy;
      const double p0 = grid.p_grid(ix, 0) + grid.p_grid(ix, 1) * cosfy +
                        grid.p_grid(ix, 2) * sinfy +
                        grid.p_grid(ix, 3) * coshy + grid.p_grid(ix, 4) * sinhy;

      // specific humidity
      const double Q = grid.q_grid(ix, 0) + grid.q_grid(ix, 1) * cosfy +
                       grid.q_grid(ix, 2) * sinfy + grid.q_grid(ix, 3) * coshy +
                       grid.q_grid(ix, 4) * sinhy;

      // lapse rate of the temperature
      g3out[k].dT = grid.dt_grid(ix, 0) + grid.dt_grid(ix, 1) * cosfy +
                    grid.dt_grid(ix, 2) * sinfy + grid.dt_grid(ix, 3) * coshy +
                    grid.dt_grid(ix, 4) * sinhy;

      // station height - grid height
      const double redh = hgt - grid.hs_grid(ix);

      // temperature at station height in Celsius
      g3out[k].T = T0 + g3out[k].dT * redh - 273.15e0;
      // printf("T = %.9f + %.9f * %.9f - %.9f\n", T0, g3out[k].dT, redh,
      // 273.15);

      // temperature lapse rate in degrees / km
      g3out[k].dT *= 1e3;

      // virtual temperature in Kelvin
      const double Tv = T0 * (1e0 + 0.6077e0 * Q);

      const double c = gm * dMtr / (Rg * Tv);

      // pressure in hPa
      g3out[k].p = (p0 * std::exp(-c * redh)) / 100e0;

      // hydrostatic and wet coefficients ah and aw
      g3out[k].ah = grid.ah_grid(ix, 0) + grid.ah_grid(ix, 1) * cosfy +
                    grid.ah_grid(ix, 2) * sinfy + grid.ah_grid(ix, 3) * coshy +
                    grid.ah_grid(ix, 4) * sinhy;
      g3out[k].aw = grid.aw_grid(ix, 0) + grid.aw_grid(ix, 1) * cosfy +
                    grid.aw_grid(ix, 2) * sinfy + grid.aw_grid(ix, 3) * coshy +
                    grid.aw_grid(ix, 4) * sinhy;

      // water vapour decrease factor la
      g3out[k].la = grid.la_grid(ix, 0) + grid.la_grid(ix, 1) * cosfy +
                    grid.la_grid(ix, 2) * sinfy + grid.la_grid(ix, 3) * coshy +
                    grid.la_grid(ix, 4) * sinhy;

      // mean temperature Tm
      g3out[k].Tm = grid.tm_grid(ix, 0) + grid.tm_grid(ix, 1) * cosfy +
                    grid.tm_grid(ix, 2) * sinfy + grid.tm_grid(ix, 3) * coshy +
                    grid.tm_grid(ix, 4) * sinhy;

      // north and east gradients (total, hydrostatic and wet)
      g3out[k].Gn_h = grid.gn_h_grid(ix, 0) + grid.gn_h_grid(ix, 1) * cosfy +
                      grid.gn_h_grid(ix, 2) * sinfy +
                      grid.gn_h_grid(ix, 3) * coshy +
                      grid.gn_h_grid(ix, 4) * sinhy;
      g3out[k].Ge_h = grid.ge_h_grid(ix, 0) + grid.ge_h_grid(ix, 1) * cosfy +
                      grid.ge_h_grid(ix, 2) * sinfy +
                      grid.ge_h_grid(ix, 3) * coshy +
                      grid.ge_h_grid(ix, 4) * sinhy;
      g3out[k].Gn_w = grid.gn_w_grid(ix, 0) + grid.gn_w_grid(ix, 1) * cosfy +
                      grid.gn_w_grid(ix, 2) * sinfy +
                      grid.gn_w_grid(ix, 3) * coshy +
                      grid.gn_w_grid(ix, 4) * sinhy;
      g3out[k].Ge_w = grid.ge_w_grid(ix, 0) + grid.ge_w_grid(ix, 1) * cosfy +
                      grid.ge_w_grid(ix, 2) * sinfy +
                      grid.ge_w_grid(ix, 3) * coshy +
                      grid.ge_w_grid(ix, 4) * sinhy;

      // water vapor pressure in hPa
      const double e0 = Q * p0 / (0.622e0 + 0.378e0 * Q) / 100e0; // on the grid
      g3out[k].e =
          e0 *
          std::pow(
              100e0 * g3out[k].p / p0,
              g3out[k].la +
                  1e0); // on the station height - (14) Askne and Nordius, 1987
    } else {
      // bilinear interpolation

      int ipod1 = ipod + sgn(diffpod);
      int ilon1 = ilon + sgn(difflon);
      if (ilon1 > ilon_max)
        ilon1 = 1;
      else if (ilon1 == 0)
        ilon1 = ilon_max;

      // get the number of the line
      indx[1] = (ipod1 - 1) * ilon_max + ilon - 1;  // along same longitude
      indx[2] = (ipod - 1) * ilon_max + ilon1 - 1;  // along same polar distance
      indx[3] = (ipod1 - 1) * ilon_max + ilon1 - 1; // diagonal
#ifdef DEBUG
      for (int i = 0; i < 4; i++)
        assert(indx[i] >= 0 && indx[i] < (int)grid.num_rows);
#endif

      // transforming ellipsoidal height to orthometric height :
      // Hortho = -N + Hell
      double undul[4] =
          {
              grid.u_grid(indx[0]),
              grid.u_grid(indx[1]),
              grid.u_grid(indx[2]),
              grid.u_grid(indx[3]),
          },
             hgt[4];
      for (int i = 0; i < 4; i++) {
        hgt[i] = /*ell[k][2]*/ hell - undul[i];
      }

      // pressure, temperature at the height of the grid
      double T0[4], p0[4];
      for (int i = 0; i < 4; i++) {
        T0[i] = grid.t_grid(indx[i], 0) + grid.t_grid(indx[i], 1) * cosfy +
                grid.t_grid(indx[i], 2) * sinfy +
                grid.t_grid(indx[i], 3) * coshy +
                grid.t_grid(indx[i], 4) * sinhy;
      }
      // printf("T = %.9f + %.9f + %.9f + %.9f \n", T0[0],T0[1],T0[2],T0[3]);
      for (int i = 0; i < 4; i++) {
        p0[i] = grid.p_grid(indx[i], 0) + grid.p_grid(indx[i], 1) * cosfy +
                grid.p_grid(indx[i], 2) * sinfy +
                grid.p_grid(indx[i], 3) * coshy +
                grid.p_grid(indx[i], 4) * sinhy;
      }
      // humidity
      double Ql[4];
      for (int i = 0; i < 4; i++) {
        Ql[i] = grid.q_grid(indx[i], 0) + grid.q_grid(indx[i], 1) * cosfy +
                grid.q_grid(indx[i], 2) * sinfy +
                grid.q_grid(indx[i], 3) * coshy +
                grid.q_grid(indx[i], 4) * sinhy;
      }

      // reduction = stationheight - gridheight
      double Hs1[4] = {grid.hs_grid(indx[0]), grid.hs_grid(indx[1]),
                       grid.hs_grid(indx[2]), grid.hs_grid(indx[3])},
             redh[4];
      for (int i = 0; i < 4; i++) {
        redh[i] = hgt[i] - Hs1[i];
      }

      // lapse rate of the temperature in degree / m
      double dTl[4];
      for (int i = 0; i < 4; i++) {
        dTl[i] = grid.dt_grid(indx[i], 0) + grid.dt_grid(indx[i], 1) * cosfy +
                 grid.dt_grid(indx[i], 2) * sinfy +
                 grid.dt_grid(indx[i], 3) * coshy +
                 grid.dt_grid(indx[i], 4) * sinhy;
      }

      // temperature reduction to station height
      double Tl[4];
      for (int i = 0; i < 4; i++) {
        Tl[i] = T0[i] + dTl[i] * redh[i] - 273.15e0;
      }

      // virtual temperature
      double Tv[4], c[4];
      for (int i = 0; i < 4; i++) {
        Tv[i] = T0[i] * (1e0 + 0.6077e0 * Ql[i]);
        c[i] = gm * dMtr / (Rg * Tv[i]);
      }

      // pressure in hPa
      double pl[4];
      for (int i = 0; i < 4; i++) {
        pl[i] = (p0[i] * std::exp(-c[i] * redh[i])) / 100e0;
      }

      // hydrostatic and wet coefficients ah and aw
      double ahl[4], awl[4];
      for (int i = 0; i < 4; i++) {
        ahl[i] = grid.ah_grid(indx[i], 0) + grid.ah_grid(indx[i], 1) * cosfy +
                 grid.ah_grid(indx[i], 2) * sinfy +
                 grid.ah_grid(indx[i], 3) * coshy +
                 grid.ah_grid(indx[i], 4) * sinhy;
      }
      for (int i = 0; i < 4; i++) {
        awl[i] = grid.aw_grid(indx[i], 0) + grid.aw_grid(indx[i], 1) * cosfy +
                 grid.aw_grid(indx[i], 2) * sinfy +
                 grid.aw_grid(indx[i], 3) * coshy +
                 grid.aw_grid(indx[i], 4) * sinhy;
      }

      // water vapour decrease factor la
      double lal[4];
      for (int i = 0; i < 4; i++) {
        lal[i] = grid.la_grid(indx[i], 0) + grid.la_grid(indx[i], 1) * cosfy +
                 grid.la_grid(indx[i], 2) * sinfy +
                 grid.la_grid(indx[i], 3) * coshy +
                 grid.la_grid(indx[i], 4) * sinhy;
      }

      // mean temperature of the water vapor Tm
      double Tml[4];
      for (int i = 0; i < 4; i++) {
        Tml[i] = grid.tm_grid(indx[i], 0) + grid.tm_grid(indx[i], 1) * cosfy +
                 grid.tm_grid(indx[i], 2) * sinfy +
                 grid.tm_grid(indx[i], 3) * coshy +
                 grid.tm_grid(indx[i], 4) * sinhy;
      }

      // north and east gradients(total, hydrostatic and wet)
      double Gn_hl[4], Ge_hl[4], Gn_wl[4], Ge_wl[4];
      for (int i = 0; i < 4; i++) {
        Gn_hl[i] = grid.gn_h_grid(indx[i], 0) +
                   grid.gn_h_grid(indx[i], 1) * cosfy +
                   grid.gn_h_grid(indx[i], 2) * sinfy +
                   grid.gn_h_grid(indx[i], 3) * coshy +
                   grid.gn_h_grid(indx[i], 4) * sinhy;
        Ge_hl[i] = grid.ge_h_grid(indx[i], 0) +
                   grid.ge_h_grid(indx[i], 1) * cosfy +
                   grid.ge_h_grid(indx[i], 2) * sinfy +
                   grid.ge_h_grid(indx[i], 3) * coshy +
                   grid.ge_h_grid(indx[i], 4) * sinhy;
        Gn_wl[i] = grid.gn_w_grid(indx[i], 0) +
                   grid.gn_w_grid(indx[i], 1) * cosfy +
                   grid.gn_w_grid(indx[i], 2) * sinfy +
                   grid.gn_w_grid(indx[i], 3) * coshy +
                   grid.gn_w_grid(indx[i], 4) * sinhy;
        Ge_wl[i] = grid.ge_w_grid(indx[i], 0) +
                   grid.ge_w_grid(indx[i], 1) * cosfy +
                   grid.ge_w_grid(indx[i], 2) * sinfy +
                   grid.ge_w_grid(indx[i], 3) * coshy +
                   grid.ge_w_grid(indx[i], 4) * sinhy;
      }

      // water vapor pressure in hPa
      double e0[4], el[4];
      for (int i = 0; i < 4; i++) {
        e0[i] =
            Ql[i] * p0[i] / (0.622e0 + 0.378e0 * Ql[i]) / 100e0; // on the grid
        el[i] =
            e0[i] *
            std::pow(
                100.e0 * pl[i] / p0[i],
                lal[i] +
                    1e0); // on the station height - (14)Askne and Nordius, 1987
      }

      double dnpod1 = std::abs(diffpod); // distance nearer podouble
      double dnpod2 = 1e0 - dnpod1;      // distance to distant podouble
      double dnlon1 = std::abs(difflon);
      double dnlon2 = 1e0 - dnlon1;

      // pressure
      double R1 = dnpod2 * pl[0] + dnpod1 * pl[1];
      double R2 = dnpod2 * pl[2] + dnpod1 * pl[3];
      g3out[k].p = dnlon2 * R1 + dnlon1 * R2;

      // temperature
      R1 = dnpod2 * Tl[0] + dnpod1 * Tl[1];
      R2 = dnpod2 * Tl[2] + dnpod1 * Tl[3];
      g3out[k].T = dnlon2 * R1 + dnlon1 * R2;

      // temperature in degree per km
      R1 = dnpod2 * dTl[0] + dnpod1 * dTl[1];
      R2 = dnpod2 * dTl[2] + dnpod1 * dTl[3];
      g3out[k].dT = (dnlon2 * R1 + dnlon1 * R2) * 1000e0;

      // water vapor pressure in hPa
      R1 = dnpod2 * el[0] + dnpod1 * el[1];
      R2 = dnpod2 * el[2] + dnpod1 * el[3];
      g3out[k].e = dnlon2 * R1 + dnlon1 * R2;

      // ah and aw
      R1 = dnpod2 * ahl[0] + dnpod1 * ahl[1];
      R2 = dnpod2 * ahl[2] + dnpod1 * ahl[3];
      g3out[k].ah = dnlon2 * R1 + dnlon1 * R2;
      R1 = dnpod2 * awl[0] + dnpod1 * awl[1];
      R2 = dnpod2 * awl[2] + dnpod1 * awl[3];
      g3out[k].aw = dnlon2 * R1 + dnlon1 * R2;

      // undulation
      R1 = dnpod2 * undul[0] + dnpod1 * undul[1];
      R2 = dnpod2 * undul[2] + dnpod1 * undul[3];
      g3out[k].undu = dnlon2 * R1 + dnlon1 * R2;

      // water vapor decrease factor la
      R1 = dnpod2 * lal[0] + dnpod1 * lal[1];
      R2 = dnpod2 * lal[2] + dnpod1 * lal[3];
      g3out[k].la = dnlon2 * R1 + dnlon1 * R2;

      // gradients
      R1 = dnpod2 * Gn_hl[0] + dnpod1 * Gn_hl[1];
      R2 = dnpod2 * Gn_hl[2] + dnpod1 * Gn_hl[3];
      g3out[k].Gn_h = dnlon2 * R1 + dnlon1 * R2;
      R1 = dnpod2 * Ge_hl[0] + dnpod1 * Ge_hl[1];
      R2 = dnpod2 * Ge_hl[2] + dnpod1 * Ge_hl[3];
      g3out[k].Ge_h = dnlon2 * R1 + dnlon1 * R2;
      R1 = dnpod2 * Gn_wl[0] + dnpod1 * Gn_wl[1];
      R2 = dnpod2 * Gn_wl[2] + dnpod1 * Gn_wl[3];
      g3out[k].Gn_w = dnlon2 * R1 + dnlon1 * R2;
      R1 = dnpod2 * Ge_wl[0] + dnpod1 * Ge_wl[1];
      R2 = dnpod2 * Ge_wl[2] + dnpod1 * Ge_wl[3];
      g3out[k].Ge_w = dnlon2 * R1 + dnlon1 * R2;

      // mean temperature of the water vapor Tm
      R1 = dnpod2 * Tml[0] + dnpod1 * Tml[1];
      R2 = dnpod2 * Tml[2] + dnpod1 * Tml[3];
      g3out[k].Tm = dnlon2 * R1 + dnlon1 * R2;
    }
  }

  return 0;
}
