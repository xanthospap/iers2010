#include "icgemio.hpp"
#include "ocean_tide.hpp"
#include <array>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace fs = std::filesystem;

namespace {

constexpr const int FES14B_MAX_DEGREE = 180;
constexpr const int FES14B_MAX_ORDER = 180;

/* (Darwin?) names of all FES2014b tidal waves; corresponding files should
 * be available to "load" the model.
 */
std::array<const char *, 34> fes14b_waves = {
    "2n2", "eps2", "j1", "k1",  "k2", "l2",  "la2", "m2",  "m3",
    "m4",  "m6",   "m8", "mf",  "mm", "mn4", "ms4", "msq", "mtm",
    "mu2", "n2",   "n4", "nu2", "o1", "om1", "om2", "p1",  "q1",
    "r2",  "s1",   "s2", "s4",  "sa", "ssa", "t2"};

int resolve_modelname(const char *modelname, const char *darwin,
                      dso::DoodsonConstituent &d) noexcept {
  /* expected string is of type: (e.g.) 'FES2014b_272.556_t2_cos' */
  const char *mdl = "FES2014b_";
  /* whitespace chars are of no importance */
  while (*modelname && *modelname == ' ')
    ++modelname;
  if (std::strncmp(mdl, modelname, 9))
    return 1;
  dso::DoodsonConstituent dnum;
  /* next part should be a valid Doodson number (in the 'IERS' format) */
  try {
    dnum = dso::resolve_iers10_doodson_string(modelname + 9);
  } catch (std::exception &) {
    return 1;
  }
  d = dnum;
  /* next part should be the Darwin string */
  const char *nam = modelname + 9 + 8; // mind the '_'
  if (std::strncmp(darwin, nam, std::strlen(darwin)))
    return 2;
  return 0;
}
} /* anonymous namespace */

dso::OceanTide dso::initFes2014bFromIcgem(const char *dir,
                                          const char *fn_generic_name,
                                          int max_degree, int max_order) {

  /* check max degree/order */
  if (max_degree < max_order || max_degree > FES14B_MAX_DEGREE ||
      max_order > FES14B_MAX_ORDER) {
    fprintf(stderr,
            "[ERROR] Invalid degree/order provided for loading FES2014b model "
            "(traceback: %s)\n",
            __func__);
    throw std::runtime_error(
        "Invalid degree/order provided for loading FES2014b model\n");
  }

  dso::OceanTide fes14b;
  fes14b.reserve(fes14b_waves.size());
  int error = 0;

  /* make sure we have a valid directory */
  fs::path gfcdir(dir);
  auto status = std::error_code{};
  if (!fs::is_directory(gfcdir, status)) {
    fprintf(stderr, "[ERROR] Failed locating directory: %s (traceback: %s)\n",
            dir, __func__);
    ++error;
  }

  if (!error) {
    char gfcfn[256];
    /* for all (major) waves of FES2014b, load two files: once with sin
     * coefficients and one with cos coeffs. Both should be located in the
     * directory provided.
     */
    for (const auto &wave_name : fes14b_waves) {
      /* resolving/adding new wave/constituent */
      dso::DoodsonConstituent doodson;
      fes14b.append_wave(doodson, max_degree, max_order);
      auto it = fes14b.waves().end() - 1;

      /* first get coefficients for sin compoment */
      try {
        std::sprintf(gfcfn, "%s/%s.%s.sin.gfc", dir, fn_generic_name,
                     wave_name);
        /* icgem instance (reads header) */
        dso::Icgem gfc(gfcfn);
        /* validate the modelname field and get the Doodson number */
        if (resolve_modelname(gfc.model_name(), wave_name, doodson)) {
          fprintf(stderr,
                  "[ERROR] Failed validating model name in gfc file %s "
                  "(traceback: %s)\n",
                  gfcfn, __func__);
          ++error;
        }
        /* parse data */
        if (!error) {
          error += gfc.parse_data(max_degree, max_order,
                                  dso::Datetime<dso::nanoseconds>::min(),
                                  it->sin_coeffs());
        }
      } catch (std::exception &e) {
        fprintf(
            stderr,
            "Failed parsing sin coefficients from file %s (traceback: %s)\n",
            gfcfn, __func__);
        ++error;
      }

      /* one more time for the cos component */
      try {
        std::sprintf(gfcfn, "%s/%s.%s.cos.gfc", dir, fn_generic_name,
                     wave_name);
        /* icgem instance (reads header) */
        dso::Icgem gfc(gfcfn);
        /* validate the modelname field and get the Doodson number */
        if (resolve_modelname(gfc.model_name(), wave_name, doodson)) {
          fprintf(stderr,
                  "[ERROR] Failed validating model name in gfc file %s "
                  "(traceback: %s)\n",
                  gfcfn, __func__);
          ++error;
        }
        /* parse data */
        if (!error) {
          error += gfc.parse_data(max_degree, max_order,
                                  dso::Datetime<dso::nanoseconds>::min(),
                                  it->cos_coeffs());
        }
      } catch (std::exception &e) {
        fprintf(
            stderr,
            "Failed parsing cos coefficients from file %s (traceback: %s)\n",
            gfcfn, __func__);
        ++error;
      }

      /* set Doodson number */
      it->doodson() = doodson;
      if (error)
        break;
    }
  }

  /* check for errors */
  if (error) {
    throw std::runtime_error("Failed loading FES14b model\n");
  }

  /* allow for computations with this degree/order */
  fes14b.resize_stokes_ceoffs(max_degree, max_order);

  return fes14b;
}
