/** @file
 *  Ocean Tide Model & computation for Geodesy.
 *
 * The easiest way to construct an Ocean Tide model, is via the published IFG 
 * files, see [2]. In this way, both major waves and admittance can be 
 * included in the model (see functions groops_ocean_atlas).
 *
 * Once a model is constructed, it can be used to compute the Stokes SH 
 * coefficients that describe the geopotential perturbation caused by the 
 * ocean tide.
 *
 * References:
 *
 * [1] Lasser, M., Meyer, U., Jäggi, A., Mayer-Gürr, T., Kvas, A., Neumayer,
 * K. H., Dahle, C., Flechtner, F., Lemoine, J.-M., Koch, I., Weigelt, M.,
 * and Flury, J.: Benchmark data for verifying background model
 * implementations in orbit and gravity field determination software,
 * Adv. Geosci., 55, 1–11, https://doi.org/10.5194/adgeo-55-1-2020, 2020.
 *
 * [2] Torsten Mayer-G ¨urr and Felix ¨Ohlinger, Gravitational effect of 
 * ocean and atmospheric tides, 
 * https://ftp.tugraz.at/outgoing/ITSG/oceanAndAtmosphericTides/README.pdf
 */

#ifndef __DSO_OCEAN_TIDE_HPP__
#define __DSO_OCEAN_TIDE_HPP__

#include "aod1b.hpp"
#include "datetime/calendar.hpp"
#include "tide_atlas.hpp"

namespace dso {

/** @class OceanTide
 *
 * A class to hold an Ocean Tide Atlas, i.e. a model to assist computing
 * the gravitational effect of oceanic tides.
 *
 * A model usually includes a series of tidal constituents (waves), each of
 * which has an individual contribution to the final effect.
 * At any given time, the gravitational effect can be represented as
 * corrections to the Stokes coefficients of a spherical harmonics expansion.
 * Any tidal constituent included in the model, should be an instance of the
 * OceanicTidalWave, so that it holds in-phase and quadrature Stokes
 * coeffcients.
 *
 * For the computation of the (corrections to the) Stokes coefficients, we
 * follow Eq. (17) of [1]
 */
class OceanTide {
private:
  /* the tidal atlas */
  TideAtlas matlas;
  /* Stokes coeffs to hold (if needed) the accumulated effect of all waves 
   * The max degree and order of this instance, should be at least equal to 
   * the maximum degree/order of the individual degree and order of the 
   * individual waves.
   */
  StokesCoeffs mcs;
  /* if we have admitance, we need to re-structure using a GroopsTideModel */
  GroopsTideModel *gmodel = nullptr;

  /** @brief Compute accumulated Stokes coefficients.
   *
   * Compute the accumulated Stokes coefficients using all available waves
   * included in the instance. The computation follows Eq. (17) of [1].
   * The computed Stokes coefficients are available in hte instance's mcs
   * member variable.
   *
   * @param[in] mjdtt Epoch for computation in TT
   * @param[in] mjdtt Epoch for computation in UT1
   * @param[in] delaunay_args An array containing Delaunay (i.e. fundamental)
   *            arguments, in the order:
   *            l  : Mean anomaly of the Moon [rad]
   *            lp : Mean anomaly of the Sun [rad]
   *            f  : L - OM [rad]
   *            d  : Mean elongation of the Moon from the Sun [rad]
   *            om : Mean longitude of the ascending node of the Moon [rad]
   * @return Always zero. Note that the computed Stokes coeffs are stored in
   *            the instance's mcs member variable.
   */
  int stokes_coeffs_impl(const MjdEpoch &mjdtt, const MjdEpoch &mjdut1,
                         const double *const delaunay_args) noexcept;
  int stokes_coeffs_with_admittance_impl(
      const MjdEpoch &mjdtt, const MjdEpoch &mjdut1,
      const double *const delaunay_args) noexcept;

public:
  /** @brief Return the underlying TideAtlas (const version) */
  const TideAtlas &atlas() const noexcept { return matlas; }

  /** @brief Return the underlying TideAtlas (non-const version) */
  TideAtlas &atlas() noexcept { return matlas; }

  /** @brief Max degree of the model (as in mcs) */
  int max_degree() const noexcept {
#ifdef DEBUG
    int maorder;
    assert(matlas.max_atlas_degree(maorder) == mcs.max_degree());
    assert(maorder == mcs.max_order());
#endif
    // return matlas.max_atlas_degree(); 
    return mcs.max_degree();
  }
  
  int max_order() const noexcept {
#ifdef DEBUG
    int maorder;
    assert(matlas.max_atlas_degree(maorder) == mcs.max_degree());
    assert(maorder == mcs.max_order());
#endif
    // return matlas.max_atlas_degree(max_order);
    return mcs.max_order();
  }

  /* @brief Return Admittance matrix (if admittance is included). */
  auto admittance_matrix() const noexcept {
    return gmodel->admittance_matrix();
  }

  auto doodson_matrix() const noexcept { return gmodel->doodson_matrix(); }
  auto num_waves() const noexcept { return gmodel->num_waves(); }
  auto num_major_waves() const noexcept { return gmodel->num_major_waves(); }

  /** @brief Constructor given an atlas and a (model) name. 
   *
   * This will construct an OcenTide model/instance, including all major 
   * waves from the Atlas, but with no admittance.
   *
   */
  OceanTide(const TideAtlas &atlas, const char *name = nullptr)
      : matlas(atlas), mcs(atlas.max_atlas_degree()) {
    int max_order;
    int max_degree = matlas.max_atlas_degree(max_order);
    if (max_order < max_degree)
      mcs.resize(max_degree, max_order);
    if (name)
      std::strcpy(matlas.name(), name);
  }

  /** @brief Constructor given an atlas, a GroopsTideModel (for admittance) and 
   * a (model) name. 
   *
   * This will construct an OcenTide model/instance, including all major waves 
   * and additionally all information to compute admittance using the 
   * GroopsTideModel passed in.
   *
   * The GroopsTideModel will be copied for this intance; hence the instance 
   * will hold its own copy and have ownership of it.
   */
  OceanTide(const TideAtlas &atlas, const GroopsTideModel &gmdl,
            const char *name = nullptr)
      : matlas(atlas), mcs(atlas.max_atlas_degree()),
        gmodel(new GroopsTideModel(gmdl)) {
    int max_order;
    int max_degree = matlas.max_atlas_degree(max_order);
    if (max_order < max_degree)
      mcs.resize(max_degree, max_order);
    if (name)
      std::strcpy(matlas.name(), name);
  }

  /** @brief Destructor; deletes GroopsTideModel if any. */
  ~OceanTide() noexcept {
    if (gmodel)
      delete gmodel;
  }

  /** @brief Copy Constructor */
  OceanTide(const OceanTide &ot) noexcept
      : matlas(ot.matlas), mcs(ot.mcs),
        gmodel((ot.gmodel) ? (new GroopsTideModel(*(ot.gmodel))) : nullptr) {}

  /** @brief Move constructor */
  OceanTide(OceanTide &&ot) noexcept
      : matlas(std::move(ot.matlas)), mcs(std::move(ot.mcs)),
        gmodel(ot.gmodel) {
    ot.gmodel = nullptr;
  }

  /** @brief Assignment operator */
  OceanTide &operator=(const OceanTide &ot) noexcept {
    if (this != &ot) {
      matlas = ot.matlas;
      mcs = ot.mcs;
      if (gmodel)
        delete gmodel;
      gmodel = (ot.gmodel) ? (new GroopsTideModel(*(ot.gmodel))) : nullptr;
    }
    return *this;
  }

  /** @brief Move assignment operator */
  OceanTide &operator=(OceanTide &&ot) noexcept {
    matlas = std::move(ot.matlas);
    mcs = std::move(ot.mcs);
    if (gmodel)
      delete gmodel;
    if (ot.gmodel) {
      gmodel = ot.gmodel;
      ot.gmodel = nullptr;
    }
    return *this;
  }

  /** @brief (Re-)Set the instance's GroopsTideModel.
   *
   * The instance will only hold a copy of the passed in GroopsTideModel.
   */
  int set_groops_admittance(const GroopsTideModel &mdl) noexcept {
    if (gmodel)
      delete gmodel;
    gmodel = new GroopsTideModel(mdl);
    return 0;
  }

  /** @brief Compute the Stokes Coefficients for a SH expansion, using the 
   * instance.
   *
   * If admittance is included (i.e. gmodel is not NULL), then except from 
   * major waves, admittance (included 'minor' waves) will also be included in 
   * the computation.
   */
  int stokes_coeffs(const MjdEpoch &mjdtt, const MjdEpoch &mjdut1,
                    const double *const delaunay_args) noexcept {
    return (gmodel) ? stokes_coeffs_with_admittance_impl(mjdtt, mjdut1,
                                                         delaunay_args)
                    : stokes_coeffs_impl(mjdtt, mjdut1, delaunay_args);
  }

  /** @brief Return the instance's Stokes Coefficients (const version) */
  const StokesCoeffs &stokes_coeffs() const noexcept { return mcs; }
  
  /** @brief Return the instance's Stokes Coefficients (non-const version) */
  StokesCoeffs &stokes_coeffs() noexcept { return mcs; }

  /** @brief Resize the instance's Stokes Coefficients (non-const version) */
  [[deprecated]]
  void resize_stokes_coeffs(int max_degree, int max_order) noexcept {
    mcs.resize(max_degree, max_order);
  }

}; /* OceanTide */

/** @brief Construct an OceanTide instance, using Groops files, see [2].
 *
 * Construct an OceanTide instance using a GROOPS files, as described in [2]. 
 * Note that this version does not include adittance.
 *
 * @param[in] fn  The <model>_001fileList.txt file of the model. For a 
 *            description of the contents/format of the file, see [2].
 * @param[in] dir The directory where the tide constituent-specific, gfc files 
 *            are located. For every major wave, two such files should 
 *            exist (for sin and cos coefficients). The filenames should 
 *            exactly match the ones specified in fn.
 * @param[in] max_degree Max Stokes Coefficients degree to parse from the 
 *            input data files. If a negative value is set, then the we are 
 *            parsing/storing up to the maximum degree set by the input files.
 * @param[in] max_order Max Stokes Coefficients order to parse from the 
 *            input data files. If a negative value is set, then the we are 
 *            parsing/storing up to the maximum order set by the input files.
 */
inline OceanTide groops_ocean_atlas(const char *fn, const char *dir,
                                    int max_degree = -1, int max_order = -1) {
  return OceanTide(groops_atlas(fn, dir, max_degree, max_order), nullptr);
}

/** @brief Construct an OceanTide instance with admittance info, using Groops 
 * files, see [2].
 *
 * Construct an OceanTide instance using a GROOPS files, as described in [2]. 
 * Note that this version does include adittance.
 *
 * @param[in] f1 The <model>_001fileList.txt file of the model. For a 
 *            description of the contents/format of the file, see [2].
 * @param[in] f2 The <model>_002doodson.txt file of the model. For a 
 *            description of the contents/format of the file, see [2].
 * @param[in] f3 The model>_003admittance.txt file of the model. For a
 *            description of the contents/format of the file, see [2].
 * @param[in] f4 The directory where the tide constituent-specific, gfc files 
 *            are located. For every major wave, two such files should 
 *            exist (for sin and cos coefficients). The filenames should 
 *            exactly match the ones specified in fn.
 * @param[in] max_degree Max Stokes Coefficients degree to parse from the 
 *            input data files. If a negative value is set, then the we are 
 *            parsing/storing up to the maximum degree set by the input files.
 * @param[in] max_order Max Stokes Coefficients order to parse from the 
 *            input data files. If a negative value is set, then the we are 
 *            parsing/storing up to the maximum order set by the input files.
 */
inline OceanTide groops_ocean_atlas(const char *f1, const char *f2,
                                    const char *f3, const char *dir,
                                    int max_degree = -1, int max_order = -1) {
  GroopsTideModel mdl;
  auto atlas{groops_atlas(f1, f2, f3, dir, mdl, max_degree, max_order)};
  return OceanTide(atlas, mdl);
}

} /* namespace dso */

#endif
