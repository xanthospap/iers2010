#ifndef __GPT3_FAST_IMPLEMENTATION_HPP__
#define __GPT3_FAST_IMPLEMENTATION_HPP__

#include "datetime/dtcalendar.hpp"
#include "geodesy/units.hpp"
#include "tropo.hpp"
#include "gpt3_parse_grid.hpp"
#include <cmath>
#ifdef USE_EXTERNAL_CONSTS
#include "gencon.hpp"
#endif
#ifdef DEBUG
#include <cstdio>
#include <cassert>
#endif

template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

namespace dso {
template<dso::gpt3::Gpt3Grid G>
int gpt3_fast(const dso::datetime<dso::nanoseconds> &t,
                     const double *lat, const double *lon, const double *hell,
                     int num_stations, int it, const gpt3::gpt3_grid<G>* gridNxN,
                     dso::gpt3_result *g3out) noexcept {
#ifdef USE_EXTERNAL_CONSTS
  constexpr double pi(DPI);
#else
  constexpr double pi(M_PI);
#endif

    // grid info
    constexpr double grid_tick = gpt3::gpt3_grid_attributes<G>::template grid_tick<double>();
    constexpr double half_grid_tick = gpt3::gpt3_grid_attributes<G>::template grid_tick<double>() / 2e0;
    constexpr int igrid_tick = gpt3::gpt3_grid_attributes<G>::template grid_tick<int>();
    printf(">> Grid tick size=%.2f, and half size=%.2f\n", grid_tick, half_grid_tick);

  // determine the GPT3 coefficients
  // mean gravity in m/s**2
  constexpr double gm = 9.80665e0;
  // molar mass of dry air in kg/mol
  constexpr double dMtr = 28.965e-3;
  // universal gas constant in J/K/mol
  constexpr double Rg = 8.3143e0;

  // get t's doy and fraction of day
  dso::ydoy_date ydoy(t.as_ydoy());
  int idoy = ydoy.__doy.as_underlying_type();
  double mjd;
  double fdoy = static_cast<double>(idoy) + std::modf(t.as_mjd(), &mjd);

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
  int bilinear = 0;

  // loop over stations
  for (int k = 0; k < num_stations; k++) {

    // only positive longitude in degrees
    double plon = dso::rad2deg(lon[k] + (lon[k] < 0) * 2e0 * pi);
    // transform to polar distance in degrees
    double ppod = dso::rad2deg(-lat[k] + pi / 2e0);

    // find the index (line in the grid file) of the nearest point
    int ipod = std::floor((ppod + igrid_tick) / igrid_tick);
    int ilon = std::floor((plon + igrid_tick) / igrid_tick);

    // normalized (to one) differences, can be positive or negative
    double diffpod = (ppod - (ipod * grid_tick - half_grid_tick)) / grid_tick;
    double difflon = (plon - (ilon * grid_tick - half_grid_tick)) / grid_tick;

    if constexpr (igrid_tick == 5) { // 5x5 grid
      if (ipod == 37)
        ipod = 36;
      if (ilon == 73)
        ilon = 1;
      else if (ilon == 0)
        ilon = 72;
    } else { // 1x1 grid
      if (ipod == 181)
        ipod = 180;
      if (ilon == 361)
        ilon = 1;
      if (ilon == 0)
        ilon = 360;
    }

    // get the number of the corresponding line
    if constexpr (igrid_tick == 5) // 5x5 grid
    indx[0] = (ipod - 1) * 72 + ilon - 1;
    else // 1x1 grid
    indx[0] = (ipod - 1) * 360 + ilon - 1;

    // near the poles: nearest neighbour interpolation, otherwise: bilinear
    bilinear = 0;
    if constexpr (igrid_tick == 5) { // 5x5 grid
      if (ppod > 2.5e0 && ppod < 177.5e0)
        bilinear = 1;
    } else { // 1x1 grid
      if (ppod > 0.5e0 && ppod < 179.5e0)
        bilinear = 1;
    }

    // case of nearest neighbourhood
    if (bilinear == 0) {

      int ix = indx[0];
#ifdef DEBUG
      assert(ix >= 0 && ix < (int)gpt3::gpt3_grid_attributes<G>::num_lines);
#endif

      // transforming ellipsoidal height to orthometric height
      g3out[k].undu = gridNxN->u_grid[ix];
      double hgt = hell[k] - g3out[k].undu;

      // pressure, temperature at the height of the grid
      double T0 = gridNxN->T_grid[ix][0] + gridNxN->T_grid[ix][1] * cosfy +
                  gridNxN->T_grid[ix][2] * sinfy + gridNxN->T_grid[ix][3] * coshy +
                  gridNxN->T_grid[ix][4] * sinhy;
      double p0 = gridNxN->p_grid[ix][0] + gridNxN->p_grid[ix][1] * cosfy +
                  gridNxN->p_grid[ix][2] * sinfy + gridNxN->p_grid[ix][3] * coshy +
                  gridNxN->p_grid[ix][4] * sinhy;

      // specific humidity
      double Q = gridNxN->Q_grid[ix][0] + gridNxN->Q_grid[ix][1] * cosfy +
                 gridNxN->Q_grid[ix][2] * sinfy + gridNxN->Q_grid[ix][3] * coshy +
                 gridNxN->Q_grid[ix][4] * sinhy;

      // lapse rate of the temperature
      g3out[k].dT = gridNxN->dT_grid[ix][0] + gridNxN->dT_grid[ix][1] * cosfy +
                    gridNxN->dT_grid[ix][2] * sinfy + gridNxN->dT_grid[ix][3] * coshy +
                    gridNxN->dT_grid[ix][4] * sinhy;

      // station height - grid height
      double redh = hgt - gridNxN->Hs_grid[ix];

      // temperature at station height in Celsius
      g3out[k].T = T0 + g3out[k].dT * redh - 273.15e0;

      // temperature lapse rate in degrees / km
      g3out[k].dT *= 1e3;

      // virtual temperature in Kelvin
      double Tv = T0 * (1e0 + 0.6077e0 * Q);

      double c = gm * dMtr / (Rg * Tv);

      // pressure in hPa
      g3out[k].p = (p0 * std::exp(-c * redh)) / 100e0;

      // hydrostatic and wet coefficients ah and aw
      g3out[k].ah = gridNxN->ah_grid[ix][0] + gridNxN->ah_grid[ix][1] * cosfy +
                    gridNxN->ah_grid[ix][2] * sinfy + gridNxN->ah_grid[ix][3] * coshy +
                    gridNxN->ah_grid[ix][4] * sinhy;
      g3out[k].aw = gridNxN->aw_grid[ix][0] + gridNxN->aw_grid[ix][1] * cosfy +
                    gridNxN->aw_grid[ix][2] * sinfy + gridNxN->aw_grid[ix][3] * coshy +
                    gridNxN->aw_grid[ix][4] * sinhy;

      // water vapour decrease factor la
      g3out[k].la = gridNxN->la_grid[ix][0] + gridNxN->la_grid[ix][1] * cosfy +
                    gridNxN->la_grid[ix][2] * sinfy + gridNxN->la_grid[ix][3] * coshy +
                    gridNxN->la_grid[ix][4] * sinhy;

      // mean temperature Tm
      g3out[k].Tm = gridNxN->Tm_grid[ix][0] + gridNxN->Tm_grid[ix][1] * cosfy +
                    gridNxN->Tm_grid[ix][2] * sinfy + gridNxN->Tm_grid[ix][3] * coshy +
                    gridNxN->Tm_grid[ix][4] * sinhy;

      // north and east gradients (total, hydrostatic and wet)
      g3out[k].Gn_h = gridNxN->Gn_h_grid[ix][0] + gridNxN->Gn_h_grid[ix][1] * cosfy +
                      gridNxN->Gn_h_grid[ix][2] * sinfy + gridNxN->Gn_h_grid[ix][3] * coshy +
                      gridNxN->Gn_h_grid[ix][4] * sinhy;
      g3out[k].Ge_h = gridNxN->Ge_h_grid[ix][0] + gridNxN->Ge_h_grid[ix][1] * cosfy +
                      gridNxN->Ge_h_grid[ix][2] * sinfy + gridNxN->Ge_h_grid[ix][3] * coshy +
                      gridNxN->Ge_h_grid[ix][4] * sinhy;
      g3out[k].Gn_w = gridNxN->Gn_w_grid[ix][0] + gridNxN->Gn_w_grid[ix][1] * cosfy +
                      gridNxN->Gn_w_grid[ix][2] * sinfy + gridNxN->Gn_w_grid[ix][3] * coshy +
                      gridNxN->Gn_w_grid[ix][4] * sinhy;
      g3out[k].Ge_w = gridNxN->Ge_w_grid[ix][0] + gridNxN->Ge_w_grid[ix][1] * cosfy +
                      gridNxN->Ge_w_grid[ix][2] * sinfy + gridNxN->Ge_w_grid[ix][3] * coshy +
                      gridNxN->Ge_w_grid[ix][4] * sinhy;

      // water vapor pressure in hPa
      double e0 = Q * p0 / (0.622e0 + 0.378e0 * Q) / 100e0; // on the grid
      g3out[k].e =
          e0 *
          std::pow(
              100e0 * g3out[k].p / p0,
              g3out[k].la +
                  1e0); // on the station height - (14) Askne and Nordius, 1987

    } else { // bilinear interpolation

      int ipod1 = ipod + sgn(diffpod);
      int ilon1 = ilon + sgn(difflon);
      if constexpr (igrid_tick == 5) { // 5x5 grid
        if (ilon1 == 73)
          ilon1 = 1;
        else if (ilon1 == 0)
          ilon1 = 72;
      } else { // 1x1 grid
        if (ilon1 == 361)
          ilon1 = 1;
        else if (ilon1 == 0)
          ilon1 = 360;
      }

      // get the number of the line
      if constexpr (igrid_tick == 5) {           // 5x5 grid
        indx[1] = (ipod1 - 1) * 72 + ilon - 1;   // along same longitude
        indx[2] = (ipod - 1) * 72 + ilon1 - 1;   // along same polar distance
        indx[3] = (ipod1 - 1) * 72 + ilon1 - 1;  // diagonal
      } else {                                   // 1x1 grid
        indx[1] = (ipod1 - 1) * 360 + ilon - 1;  // along same longitude
        indx[2] = (ipod - 1) * 360 + ilon1 - 1;  // along same polar distance
        indx[3] = (ipod1 - 1) * 360 + ilon1 - 1; // diagonal
      }
#ifdef DEBUG
      for (int i = 0; i < 4; i++)
        assert(indx[i] >= 0 && indx[i] < (int)gpt3::gpt3_grid_attributes<G>::num_lines);
#endif

      // transforming ellipsoidal height to orthometric height :
      // Hortho = -N + Hell
      double undul[4] =
          {
              gridNxN->u_grid[indx[0]],
              gridNxN->u_grid[indx[1]],
              gridNxN->u_grid[indx[2]],
              gridNxN->u_grid[indx[3]],
          },
             hgt[4];
      for (int i = 0; i < 4; i++) {
        hgt[i] = hell[k] - undul[i];
      }

      // pressure, temperature at the height of the grid
      double T0[4], p0[4];
      for (int i = 0; i < 4; i++) {
        T0[i] = gridNxN->T_grid[indx[i]][0] + gridNxN->T_grid[indx[i]][1] * cosfy +
                gridNxN->T_grid[indx[i]][2] * sinfy + gridNxN->T_grid[indx[i]][3] * coshy +
                gridNxN->T_grid[indx[i]][4] * sinhy;
      }
      for (int i = 0; i < 4; i++) {
        p0[i] = gridNxN->p_grid[indx[i]][0] + gridNxN->p_grid[indx[i]][1] * cosfy +
                gridNxN->p_grid[indx[i]][2] * sinfy + gridNxN->p_grid[indx[i]][3] * coshy +
                gridNxN->p_grid[indx[i]][4] * sinhy;
      }
      // humidity
      double Ql[4];
      for (int i = 0; i < 4; i++) {
        Ql[i] = gridNxN->Q_grid[indx[i]][0] + gridNxN->Q_grid[indx[i]][1] * cosfy +
                gridNxN->Q_grid[indx[i]][2] * sinfy + gridNxN->Q_grid[indx[i]][3] * coshy +
                gridNxN->Q_grid[indx[i]][4] * sinhy;
      }

      // reduction = stationheight - gridheight
      double Hs1[4] = {gridNxN->Hs_grid[indx[0]], gridNxN->Hs_grid[indx[1]],
                       gridNxN->Hs_grid[indx[2]], gridNxN->Hs_grid[indx[3]]},
             redh[4];
      for (int i = 0; i < 4; i++) {
        redh[i] = hgt[i] - Hs1[i];
      }

      // lapse rate of the temperature in degree / m
      double dTl[4];
      for (int i = 0; i < 4; i++) {
        dTl[i] = gridNxN->dT_grid[indx[i]][0] + gridNxN->dT_grid[indx[i]][1] * cosfy +
                 gridNxN->dT_grid[indx[i]][2] * sinfy + gridNxN->dT_grid[indx[i]][3] * coshy +
                 gridNxN->dT_grid[indx[i]][4] * sinhy;
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
        ahl[i] = gridNxN->ah_grid[indx[i]][0] + gridNxN->ah_grid[indx[i]][1] * cosfy +
                 gridNxN->ah_grid[indx[i]][2] * sinfy + gridNxN->ah_grid[indx[i]][3] * coshy +
                 gridNxN->ah_grid[indx[i]][4] * sinhy;
      }
      for (int i = 0; i < 4; i++) {
        awl[i] = gridNxN->aw_grid[indx[i]][0] + gridNxN->aw_grid[indx[i]][1] * cosfy +
                 gridNxN->aw_grid[indx[i]][2] * sinfy + gridNxN->aw_grid[indx[i]][3] * coshy +
                 gridNxN->aw_grid[indx[i]][4] * sinhy;
      }

      // water vapour decrease factor la
      double lal[4];
      for (int i = 0; i < 4; i++) {
        lal[i] = gridNxN->la_grid[indx[i]][0] + gridNxN->la_grid[indx[i]][1] * cosfy +
                 gridNxN->la_grid[indx[i]][2] * sinfy + gridNxN->la_grid[indx[i]][3] * coshy +
                 gridNxN->la_grid[indx[i]][4] * sinhy;
      }

      // mean temperature of the water vapor Tm
      double Tml[4];
      for (int i = 0; i < 4; i++) {
        Tml[i] = gridNxN->Tm_grid[indx[i]][0] + gridNxN->Tm_grid[indx[i]][1] * cosfy +
                 gridNxN->Tm_grid[indx[i]][2] * sinfy + gridNxN->Tm_grid[indx[i]][3] * coshy +
                 gridNxN->Tm_grid[indx[i]][4] * sinhy;
      }

      // north and east gradients(total, hydrostatic and wet)
      double Gn_hl[4], Ge_hl[4], Gn_wl[4], Ge_wl[4];
      for (int i = 0; i < 4; i++) {
        Gn_hl[i] = gridNxN->Gn_h_grid[indx[i]][0] + gridNxN->Gn_h_grid[indx[i]][1] * cosfy +
                   gridNxN->Gn_h_grid[indx[i]][2] * sinfy +
                   gridNxN->Gn_h_grid[indx[i]][3] * coshy +
                   gridNxN->Gn_h_grid[indx[i]][4] * sinhy;
        Ge_hl[i] = gridNxN->Ge_h_grid[indx[i]][0] + gridNxN->Ge_h_grid[indx[i]][1] * cosfy +
                   gridNxN->Ge_h_grid[indx[i]][2] * sinfy +
                   gridNxN->Ge_h_grid[indx[i]][3] * coshy +
                   gridNxN->Ge_h_grid[indx[i]][4] * sinhy;
        Gn_wl[i] = gridNxN->Gn_w_grid[indx[i]][0] + gridNxN->Gn_w_grid[indx[i]][1] * cosfy +
                   gridNxN->Gn_w_grid[indx[i]][2] * sinfy +
                   gridNxN->Gn_w_grid[indx[i]][3] * coshy +
                   gridNxN->Gn_w_grid[indx[i]][4] * sinhy;
        Ge_wl[i] = gridNxN->Ge_w_grid[indx[i]][0] + gridNxN->Ge_w_grid[indx[i]][1] * cosfy +
                   gridNxN->Ge_w_grid[indx[i]][2] * sinfy +
                   gridNxN->Ge_w_grid[indx[i]][3] * coshy +
                   gridNxN->Ge_w_grid[indx[i]][4] * sinhy;
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

template<dso::gpt3::Gpt3Grid G>
int gpt3_fast(const dso::datetime<dso::nanoseconds> &t,
                     const double *lat, const double *lon, const double *hell,
                     int num_stations, int it, const char *grid_file,
                     dso::gpt3_result *g3out) noexcept {
  printf(">> called gpt3_fast(#1), with grid file=%s\n", grid_file);
  gpt3::gpt3_grid<G> gridNxN;
  printf(">> Size of grid instance is: %lu\n", sizeof(gridNxN));

  if (gpt3::parse_gpt3_grid<G>(grid_file, &gridNxN)) {
    fprintf(stderr,
            "[ERROR] Failed parsing gpt3 grid file! (traceback: %s)\n",
            __func__);
    return 15;
  }

  return gpt3_fast(t, lat, lon, hell, num_stations, it, &gridNxN, g3out);
}

}// dso
#endif
