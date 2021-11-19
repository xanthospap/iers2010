#include "ggdatetime/dtcalendar.hpp"
#include "ggeodesy/units.hpp"
#include "tropo.hpp"
#include <cmath>
#ifdef USE_EXTERNAL_CONSTS
#include "gencon.hpp"
#endif
#ifdef DEBUG
#include <cassert>
#endif

template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

int dso::gpt3_5_fast(const dso::datetime<dso::nanoseconds> &t,
                     const double *lat, const double *lon, const double *hell,
                     int num_stations, int it, const char *grid_file,
                     dso::gpt3_result *g3out) noexcept {
/*
#ifdef USE_EXTERNAL_CONSTS
  constexpr double pi(DPI);
#else
  constexpr double pi(M_PI);
#endif

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
*/
  gpt3::gpt3_5_grid grid5x5;
  if (gpt3::parse_gpt3_5_grid(grid_file, &grid5x5)) {
    fprintf(stderr,
            "[ERROR] Failed parsing gpt3_5 grid file! (traceback: %s)\n",
            __func__);
    return 15;
  }

  return gpt3_5_fast(t, lat, lon, hell, num_stations, it, &grid5x5, g3out);
/*
  // loop over stations
  for (int k = 0; k < num_stations; k++) {

    // only positive longitude in degrees
    double plon = dso::rad2deg(lon[k] + (lon[k] < 0) * 2e0 * pi);
    // transform to polar distance in degrees
    double ppod = dso::rad2deg(-lat[k] + pi / 2e0);

    // find the index (line in the grid file) of the nearest point
    int ipod = std::floor((ppod + 5) / 5);
    int ilon = std::floor((plon + 5) / 5);

    // normalized (to one) differences, can be positive or negative
    double diffpod = (ppod - (ipod * 5e0 - 2.5e0)) / 5e0;
    double difflon = (plon - (ilon * 5e0 - 2.5e0)) / 5e0;
    if (ipod == 37)
      ipod = 36;
    if (ilon == 73)
      ilon = 1;
    else if (ilon == 0)
      ilon = 72;

    // get the number of the corresponding line
    indx[0] = (ipod - 1) * 72 + ilon - 1;

    // near the poles: nearest neighbour interpolation, otherwise: bilinear
    bilinear = 0;
    if (ppod > 2.5 && ppod < 177.5)
      bilinear = 1;

    // case of nearest neighbourhood
    if (bilinear == 0) {

      int ix = indx[0];
#ifdef DEBUG
      assert(ix >= 0 && ix < gpt3::GPT3_5_GRID_LINES);
#endif

      // transforming ellipsoidal height to orthometric height
      g3out[k].undu = g.u_grid[ix];
      double hgt = hell[k] - g3out[k].undu;

      // pressure, temperature at the height of the grid
      double T0 = g.T_grid[ix][0] + g.T_grid[ix][1] * cosfy +
                  g.T_grid[ix][2] * sinfy + g.T_grid[ix][3] * coshy +
                  g.T_grid[ix][4] * sinhy;
      double p0 = g.p_grid[ix][0] + g.p_grid[ix][1] * cosfy +
                  g.p_grid[ix][2] * sinfy + g.p_grid[ix][3] * coshy +
                  g.p_grid[ix][4] * sinhy;

      // specific humidity
      double Q = g.Q_grid[ix][0] + g.Q_grid[ix][1] * cosfy +
                 g.Q_grid[ix][2] * sinfy + g.Q_grid[ix][3] * coshy +
                 g.Q_grid[ix][4] * sinhy;

      // lapse rate of the temperature
      g3out[k].dT = g.dT_grid[ix][0] + g.dT_grid[ix][1] * cosfy +
                    g.dT_grid[ix][2] * sinfy + g.dT_grid[ix][3] * coshy +
                    g.dT_grid[ix][4] * sinhy;

      // station height - grid height
      double redh = hgt - g.Hs_grid[ix];

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
      g3out[k].ah = g.ah_grid[ix][0] + g.ah_grid[ix][1] * cosfy +
                    g.ah_grid[ix][2] * sinfy + g.ah_grid[ix][3] * coshy +
                    g.ah_grid[ix][4] * sinhy;
      g3out[k].aw = g.aw_grid[ix][0] + g.aw_grid[ix][1] * cosfy +
                    g.aw_grid[ix][2] * sinfy + g.aw_grid[ix][3] * coshy +
                    g.aw_grid[ix][4] * sinhy;

      // water vapour decrease factor la
      g3out[k].la = g.la_grid[ix][0] + g.la_grid[ix][1] * cosfy +
                    g.la_grid[ix][2] * sinfy + g.la_grid[ix][3] * coshy +
                    g.la_grid[ix][4] * sinhy;

      // mean temperature Tm
      g3out[k].Tm = g.Tm_grid[ix][0] + g.Tm_grid[ix][1] * cosfy +
                    g.Tm_grid[ix][2] * sinfy + g.Tm_grid[ix][3] * coshy +
                    g.Tm_grid[ix][4] * sinhy;

      // north and east gradients (total, hydrostatic and wet)
      g3out[k].Gn_h = g.Gn_h_grid[ix][0] + g.Gn_h_grid[ix][1] * cosfy +
                      g.Gn_h_grid[ix][2] * sinfy + g.Gn_h_grid[ix][3] * coshy +
                      g.Gn_h_grid[ix][4] * sinhy;
      g3out[k].Ge_h = g.Ge_h_grid[ix][0] + g.Ge_h_grid[ix][1] * cosfy +
                      g.Ge_h_grid[ix][2] * sinfy + g.Ge_h_grid[ix][3] * coshy +
                      g.Ge_h_grid[ix][4] * sinhy;
      g3out[k].Gn_w = g.Gn_w_grid[ix][0] + g.Gn_w_grid[ix][1] * cosfy +
                      g.Gn_w_grid[ix][2] * sinfy + g.Gn_w_grid[ix][3] * coshy +
                      g.Gn_w_grid[ix][4] * sinhy;
      g3out[k].Ge_w = g.Ge_w_grid[ix][0] + g.Ge_w_grid[ix][1] * cosfy +
                      g.Ge_w_grid[ix][2] * sinfy + g.Ge_w_grid[ix][3] * coshy +
                      g.Ge_w_grid[ix][4] * sinhy;

      // water vapor pressure in hPa
      double e0 = Q * p0 / (0.622e0 + 0.378e0 * Q) / 100e0; // on the grid
      g3out[k].e =
          e0 *
          std::pow(
              100e0 * g3out[k].p / p0,
              g3out[k].la +
                  1e0); // on the station height - (14) Askne and Nordius, 1987

    } else { //% bilinear interpolation

      int ipod1 = ipod + sgn(diffpod);
      int ilon1 = ilon + sgn(difflon);
      if (ilon1 == 73)
        ilon1 = 1;
      else if (ilon1 == 0)
        ilon1 = 72;

      // get the number of the line
      indx[1] = (ipod1 - 1) * 72 + ilon - 1;  // along same longitude
      indx[2] = (ipod - 1) * 72 + ilon1 - 1;  // along same polar distance
      indx[3] = (ipod1 - 1) * 72 + ilon1 - 1; // diagonal
#ifdef DEBUG
      for (int i = 0; i < 4; i++)
        assert(indx[i] >= 0 && indx[i] < gpt3::GPT3_5_GRID_LINES);
#endif

      // transforming ellipsoidal height to orthometric height :
      // Hortho = -N + Hell
      double undul[4] =
          {
              g.u_grid[indx[0]],
              g.u_grid[indx[1]],
              g.u_grid[indx[2]],
              g.u_grid[indx[3]],
          },
             hgt[4];
      for (int i = 0; i < 4; i++) {
        hgt[i] = hell[k] - undul[i];
      }

      // pressure, temperature at the height of the grid
      double T0[4], p0[4];
      for (int i = 0; i < 4; i++) {
        T0[i] = g.T_grid[indx[i]][0] + g.T_grid[indx[i]][1] * cosfy +
                g.T_grid[indx[i]][2] * sinfy + g.T_grid[indx[i]][3] * coshy +
                g.T_grid[indx[i]][4] * sinhy;
      }
      for (int i = 0; i < 4; i++) {
        p0[i] = g.p_grid[indx[i]][0] + g.p_grid[indx[i]][1] * cosfy +
                g.p_grid[indx[i]][2] * sinfy + g.p_grid[indx[i]][3] * coshy +
                g.p_grid[indx[i]][4] * sinhy;
      }
      // humidity
      double Ql[4];
      for (int i = 0; i < 4; i++) {
        Ql[i] = g.Q_grid[indx[i]][0] + g.Q_grid[indx[i]][1] * cosfy +
                g.Q_grid[indx[i]][2] * sinfy + g.Q_grid[indx[i]][3] * coshy +
                g.Q_grid[indx[i]][4] * sinhy;
      }

      // reduction = stationheight - gridheight
      double Hs1[4] = {g.Hs_grid[indx[0]], g.Hs_grid[indx[1]],
                       g.Hs_grid[indx[2]], g.Hs_grid[indx[3]]},
             redh[4];
      for (int i = 0; i < 4; i++) {
        redh[i] = hgt[i] - Hs1[i];
      }

      // lapse rate of the temperature in degree / m
      double dTl[4];
      for (int i = 0; i < 4; i++) {
        dTl[i] = g.dT_grid[indx[i]][0] + g.dT_grid[indx[i]][1] * cosfy +
                 g.dT_grid[indx[i]][2] * sinfy + g.dT_grid[indx[i]][3] * coshy +
                 g.dT_grid[indx[i]][4] * sinhy;
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
        ahl[i] = g.ah_grid[indx[i]][0] + g.ah_grid[indx[i]][1] * cosfy +
                 g.ah_grid[indx[i]][2] * sinfy + g.ah_grid[indx[i]][3] * coshy +
                 g.ah_grid[indx[i]][4] * sinhy;
      }
      for (int i = 0; i < 4; i++) {
        awl[i] = g.aw_grid[indx[i]][0] + g.aw_grid[indx[i]][1] * cosfy +
                 g.aw_grid[indx[i]][2] * sinfy + g.aw_grid[indx[i]][3] * coshy +
                 g.aw_grid[indx[i]][4] * sinhy;
      }

      // water vapour decrease factor la
      double lal[4];
      for (int i = 0; i < 4; i++) {
        lal[i] = g.la_grid[indx[i]][0] + g.la_grid[indx[i]][1] * cosfy +
                 g.la_grid[indx[i]][2] * sinfy + g.la_grid[indx[i]][3] * coshy +
                 g.la_grid[indx[i]][4] * sinhy;
      }

      // mean temperature of the water vapor Tm
      double Tml[4];
      for (int i = 0; i < 4; i++) {
        Tml[i] = g.Tm_grid[indx[i]][0] + g.Tm_grid[indx[i]][1] * cosfy +
                 g.Tm_grid[indx[i]][2] * sinfy + g.Tm_grid[indx[i]][3] * coshy +
                 g.Tm_grid[indx[i]][4] * sinhy;
      }

      // north and east gradients(total, hydrostatic and wet)
      double Gn_hl[4], Ge_hl[4], Gn_wl[4], Ge_wl[4];
      for (int i = 0; i < 4; i++) {
        Gn_hl[i] = g.Gn_h_grid[indx[i]][0] + g.Gn_h_grid[indx[i]][1] * cosfy +
                   g.Gn_h_grid[indx[i]][2] * sinfy +
                   g.Gn_h_grid[indx[i]][3] * coshy +
                   g.Gn_h_grid[indx[i]][4] * sinhy;
        Ge_hl[i] = g.Ge_h_grid[indx[i]][0] + g.Ge_h_grid[indx[i]][1] * cosfy +
                   g.Ge_h_grid[indx[i]][2] * sinfy +
                   g.Ge_h_grid[indx[i]][3] * coshy +
                   g.Ge_h_grid[indx[i]][4] * sinhy;
        Gn_wl[i] = g.Gn_w_grid[indx[i]][0] + g.Gn_w_grid[indx[i]][1] * cosfy +
                   g.Gn_w_grid[indx[i]][2] * sinfy +
                   g.Gn_w_grid[indx[i]][3] * coshy +
                   g.Gn_w_grid[indx[i]][4] * sinhy;
        Ge_wl[i] = g.Ge_w_grid[indx[i]][0] + g.Ge_w_grid[indx[i]][1] * cosfy +
                   g.Ge_w_grid[indx[i]][2] * sinfy +
                   g.Ge_w_grid[indx[i]][3] * coshy +
                   g.Ge_w_grid[indx[i]][4] * sinhy;
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
  */
}

int dso::gpt3_5_fast(const dso::datetime<dso::nanoseconds> &t,
                     const double *lat, const double *lon, const double *hell,
                     int num_stations, int it, const gpt3::gpt3_5_grid* grid5x5,
                     dso::gpt3_result *g3out) noexcept {
#ifdef USE_EXTERNAL_CONSTS
  constexpr double pi(DPI);
#else
  constexpr double pi(M_PI);
#endif

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

  /*gpt3::gpt3_5_grid g;
  if (gpt3::parse_gpt3_5_grid(grid_file, &g)) {
    fprintf(stderr,
            "[ERROR] Failed parsing gpt3_5 grid file! (traceback: %s)\n",
            __func__);
    return 15;
  }*/

  // loop over stations
  for (int k = 0; k < num_stations; k++) {

    // only positive longitude in degrees
    double plon = dso::rad2deg(lon[k] + (lon[k] < 0) * 2e0 * pi);
    // transform to polar distance in degrees
    double ppod = dso::rad2deg(-lat[k] + pi / 2e0);

    // find the index (line in the grid file) of the nearest point
    int ipod = std::floor((ppod + 5) / 5);
    int ilon = std::floor((plon + 5) / 5);

    // normalized (to one) differences, can be positive or negative
    double diffpod = (ppod - (ipod * 5e0 - 2.5e0)) / 5e0;
    double difflon = (plon - (ilon * 5e0 - 2.5e0)) / 5e0;
    if (ipod == 37)
      ipod = 36;
    if (ilon == 73)
      ilon = 1;
    else if (ilon == 0)
      ilon = 72;

    // get the number of the corresponding line
    indx[0] = (ipod - 1) * 72 + ilon - 1;

    // near the poles: nearest neighbour interpolation, otherwise: bilinear
    bilinear = 0;
    if (ppod > 2.5 && ppod < 177.5)
      bilinear = 1;

    // case of nearest neighbourhood
    if (bilinear == 0) {

      int ix = indx[0];
#ifdef DEBUG
      assert(ix >= 0 && ix < gpt3::GPT3_5_GRID_LINES);
#endif

      // transforming ellipsoidal height to orthometric height
      g3out[k].undu = grid5x5->u_grid[ix];
      double hgt = hell[k] - g3out[k].undu;

      // pressure, temperature at the height of the grid
      double T0 = grid5x5->T_grid[ix][0] + grid5x5->T_grid[ix][1] * cosfy +
                  grid5x5->T_grid[ix][2] * sinfy + grid5x5->T_grid[ix][3] * coshy +
                  grid5x5->T_grid[ix][4] * sinhy;
      double p0 = grid5x5->p_grid[ix][0] + grid5x5->p_grid[ix][1] * cosfy +
                  grid5x5->p_grid[ix][2] * sinfy + grid5x5->p_grid[ix][3] * coshy +
                  grid5x5->p_grid[ix][4] * sinhy;

      // specific humidity
      double Q = grid5x5->Q_grid[ix][0] + grid5x5->Q_grid[ix][1] * cosfy +
                 grid5x5->Q_grid[ix][2] * sinfy + grid5x5->Q_grid[ix][3] * coshy +
                 grid5x5->Q_grid[ix][4] * sinhy;

      // lapse rate of the temperature
      g3out[k].dT = grid5x5->dT_grid[ix][0] + grid5x5->dT_grid[ix][1] * cosfy +
                    grid5x5->dT_grid[ix][2] * sinfy + grid5x5->dT_grid[ix][3] * coshy +
                    grid5x5->dT_grid[ix][4] * sinhy;

      // station height - grid height
      double redh = hgt - grid5x5->Hs_grid[ix];

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
      g3out[k].ah = grid5x5->ah_grid[ix][0] + grid5x5->ah_grid[ix][1] * cosfy +
                    grid5x5->ah_grid[ix][2] * sinfy + grid5x5->ah_grid[ix][3] * coshy +
                    grid5x5->ah_grid[ix][4] * sinhy;
      g3out[k].aw = grid5x5->aw_grid[ix][0] + grid5x5->aw_grid[ix][1] * cosfy +
                    grid5x5->aw_grid[ix][2] * sinfy + grid5x5->aw_grid[ix][3] * coshy +
                    grid5x5->aw_grid[ix][4] * sinhy;

      // water vapour decrease factor la
      g3out[k].la = grid5x5->la_grid[ix][0] + grid5x5->la_grid[ix][1] * cosfy +
                    grid5x5->la_grid[ix][2] * sinfy + grid5x5->la_grid[ix][3] * coshy +
                    grid5x5->la_grid[ix][4] * sinhy;

      // mean temperature Tm
      g3out[k].Tm = grid5x5->Tm_grid[ix][0] + grid5x5->Tm_grid[ix][1] * cosfy +
                    grid5x5->Tm_grid[ix][2] * sinfy + grid5x5->Tm_grid[ix][3] * coshy +
                    grid5x5->Tm_grid[ix][4] * sinhy;

      // north and east gradients (total, hydrostatic and wet)
      g3out[k].Gn_h = grid5x5->Gn_h_grid[ix][0] + grid5x5->Gn_h_grid[ix][1] * cosfy +
                      grid5x5->Gn_h_grid[ix][2] * sinfy + grid5x5->Gn_h_grid[ix][3] * coshy +
                      grid5x5->Gn_h_grid[ix][4] * sinhy;
      g3out[k].Ge_h = grid5x5->Ge_h_grid[ix][0] + grid5x5->Ge_h_grid[ix][1] * cosfy +
                      grid5x5->Ge_h_grid[ix][2] * sinfy + grid5x5->Ge_h_grid[ix][3] * coshy +
                      grid5x5->Ge_h_grid[ix][4] * sinhy;
      g3out[k].Gn_w = grid5x5->Gn_w_grid[ix][0] + grid5x5->Gn_w_grid[ix][1] * cosfy +
                      grid5x5->Gn_w_grid[ix][2] * sinfy + grid5x5->Gn_w_grid[ix][3] * coshy +
                      grid5x5->Gn_w_grid[ix][4] * sinhy;
      g3out[k].Ge_w = grid5x5->Ge_w_grid[ix][0] + grid5x5->Ge_w_grid[ix][1] * cosfy +
                      grid5x5->Ge_w_grid[ix][2] * sinfy + grid5x5->Ge_w_grid[ix][3] * coshy +
                      grid5x5->Ge_w_grid[ix][4] * sinhy;

      // water vapor pressure in hPa
      double e0 = Q * p0 / (0.622e0 + 0.378e0 * Q) / 100e0; // on the grid
      g3out[k].e =
          e0 *
          std::pow(
              100e0 * g3out[k].p / p0,
              g3out[k].la +
                  1e0); // on the station height - (14) Askne and Nordius, 1987

    } else { //% bilinear interpolation

      int ipod1 = ipod + sgn(diffpod);
      int ilon1 = ilon + sgn(difflon);
      if (ilon1 == 73)
        ilon1 = 1;
      else if (ilon1 == 0)
        ilon1 = 72;

      // get the number of the line
      indx[1] = (ipod1 - 1) * 72 + ilon - 1;  // along same longitude
      indx[2] = (ipod - 1) * 72 + ilon1 - 1;  // along same polar distance
      indx[3] = (ipod1 - 1) * 72 + ilon1 - 1; // diagonal
#ifdef DEBUG
      for (int i = 0; i < 4; i++)
        assert(indx[i] >= 0 && indx[i] < gpt3::GPT3_5_GRID_LINES);
#endif

      // transforming ellipsoidal height to orthometric height :
      // Hortho = -N + Hell
      double undul[4] =
          {
              grid5x5->u_grid[indx[0]],
              grid5x5->u_grid[indx[1]],
              grid5x5->u_grid[indx[2]],
              grid5x5->u_grid[indx[3]],
          },
             hgt[4];
      for (int i = 0; i < 4; i++) {
        hgt[i] = hell[k] - undul[i];
      }

      // pressure, temperature at the height of the grid
      double T0[4], p0[4];
      for (int i = 0; i < 4; i++) {
        T0[i] = grid5x5->T_grid[indx[i]][0] + grid5x5->T_grid[indx[i]][1] * cosfy +
                grid5x5->T_grid[indx[i]][2] * sinfy + grid5x5->T_grid[indx[i]][3] * coshy +
                grid5x5->T_grid[indx[i]][4] * sinhy;
      }
      for (int i = 0; i < 4; i++) {
        p0[i] = grid5x5->p_grid[indx[i]][0] + grid5x5->p_grid[indx[i]][1] * cosfy +
                grid5x5->p_grid[indx[i]][2] * sinfy + grid5x5->p_grid[indx[i]][3] * coshy +
                grid5x5->p_grid[indx[i]][4] * sinhy;
      }
      // humidity
      double Ql[4];
      for (int i = 0; i < 4; i++) {
        Ql[i] = grid5x5->Q_grid[indx[i]][0] + grid5x5->Q_grid[indx[i]][1] * cosfy +
                grid5x5->Q_grid[indx[i]][2] * sinfy + grid5x5->Q_grid[indx[i]][3] * coshy +
                grid5x5->Q_grid[indx[i]][4] * sinhy;
      }

      // reduction = stationheight - gridheight
      double Hs1[4] = {grid5x5->Hs_grid[indx[0]], grid5x5->Hs_grid[indx[1]],
                       grid5x5->Hs_grid[indx[2]], grid5x5->Hs_grid[indx[3]]},
             redh[4];
      for (int i = 0; i < 4; i++) {
        redh[i] = hgt[i] - Hs1[i];
      }

      // lapse rate of the temperature in degree / m
      double dTl[4];
      for (int i = 0; i < 4; i++) {
        dTl[i] = grid5x5->dT_grid[indx[i]][0] + grid5x5->dT_grid[indx[i]][1] * cosfy +
                 grid5x5->dT_grid[indx[i]][2] * sinfy + grid5x5->dT_grid[indx[i]][3] * coshy +
                 grid5x5->dT_grid[indx[i]][4] * sinhy;
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
        ahl[i] = grid5x5->ah_grid[indx[i]][0] + grid5x5->ah_grid[indx[i]][1] * cosfy +
                 grid5x5->ah_grid[indx[i]][2] * sinfy + grid5x5->ah_grid[indx[i]][3] * coshy +
                 grid5x5->ah_grid[indx[i]][4] * sinhy;
      }
      for (int i = 0; i < 4; i++) {
        awl[i] = grid5x5->aw_grid[indx[i]][0] + grid5x5->aw_grid[indx[i]][1] * cosfy +
                 grid5x5->aw_grid[indx[i]][2] * sinfy + grid5x5->aw_grid[indx[i]][3] * coshy +
                 grid5x5->aw_grid[indx[i]][4] * sinhy;
      }

      // water vapour decrease factor la
      double lal[4];
      for (int i = 0; i < 4; i++) {
        lal[i] = grid5x5->la_grid[indx[i]][0] + grid5x5->la_grid[indx[i]][1] * cosfy +
                 grid5x5->la_grid[indx[i]][2] * sinfy + grid5x5->la_grid[indx[i]][3] * coshy +
                 grid5x5->la_grid[indx[i]][4] * sinhy;
      }

      // mean temperature of the water vapor Tm
      double Tml[4];
      for (int i = 0; i < 4; i++) {
        Tml[i] = grid5x5->Tm_grid[indx[i]][0] + grid5x5->Tm_grid[indx[i]][1] * cosfy +
                 grid5x5->Tm_grid[indx[i]][2] * sinfy + grid5x5->Tm_grid[indx[i]][3] * coshy +
                 grid5x5->Tm_grid[indx[i]][4] * sinhy;
      }

      // north and east gradients(total, hydrostatic and wet)
      double Gn_hl[4], Ge_hl[4], Gn_wl[4], Ge_wl[4];
      for (int i = 0; i < 4; i++) {
        Gn_hl[i] = grid5x5->Gn_h_grid[indx[i]][0] + grid5x5->Gn_h_grid[indx[i]][1] * cosfy +
                   grid5x5->Gn_h_grid[indx[i]][2] * sinfy +
                   grid5x5->Gn_h_grid[indx[i]][3] * coshy +
                   grid5x5->Gn_h_grid[indx[i]][4] * sinhy;
        Ge_hl[i] = grid5x5->Ge_h_grid[indx[i]][0] + grid5x5->Ge_h_grid[indx[i]][1] * cosfy +
                   grid5x5->Ge_h_grid[indx[i]][2] * sinfy +
                   grid5x5->Ge_h_grid[indx[i]][3] * coshy +
                   grid5x5->Ge_h_grid[indx[i]][4] * sinhy;
        Gn_wl[i] = grid5x5->Gn_w_grid[indx[i]][0] + grid5x5->Gn_w_grid[indx[i]][1] * cosfy +
                   grid5x5->Gn_w_grid[indx[i]][2] * sinfy +
                   grid5x5->Gn_w_grid[indx[i]][3] * coshy +
                   grid5x5->Gn_w_grid[indx[i]][4] * sinhy;
        Ge_wl[i] = grid5x5->Ge_w_grid[indx[i]][0] + grid5x5->Ge_w_grid[indx[i]][1] * cosfy +
                   grid5x5->Ge_w_grid[indx[i]][2] * sinfy +
                   grid5x5->Ge_w_grid[indx[i]][3] * coshy +
                   grid5x5->Ge_w_grid[indx[i]][4] * sinhy;
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
