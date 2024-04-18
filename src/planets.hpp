/** @file
 *
 * Fundamental operations for the determination of planetary states (i.e.
 * position & velocity) for major planets (and the Moon) of the solar system.
 *
 * This file also includes an API to the JPL DE ephemeris, via the CSPICE
 * Toolkit.
 */

#ifndef __DSO_PLANETARY_STATE_AND_JPLDE_HPP__
#define __DSO_PLANETARY_STATE_AND_JPLDE_HPP__

#include "cppspice/SpiceUsr.h"
#include "datetime/calendar.hpp"
#include "eigen3/Eigen/Eigen"

namespace dso {

/* Enum class to denote kinda planets, for easy interaction with CSPICE */
enum class Planet : char {
  EARTH,
  MOON,
  SUN,
  MERCURY,
  VENUS,
  MARS,
  JUPITER,
  SATURN
};

namespace cspice {

/* @brief Transform a Planet enum to its NAIF Id.
 *
 * @param[in] p A planet included in Planet enum class
 * @param[out] id The planet's NAIF ID (to be used in the library's API)
 * @return If a value other than zero is returned, then we couldn't match 
 *         an (int) id to the given planet.
 */
int planet_to_naif_id(dso::Planet p, int &id) noexcept;

/* @brief Transform a Planet enum to its NAIF Id and retuen it as string.
 *
 * @param[in]  p A planet included in Planet enum class
 * @param[out] buf The planet's NAIF ID (to be used in the library's API), 
 *             returned as string. E.g. Earth's NAIF id is 399; hence, at 
 *             return, buf will hold the string '399\0' (the returned string 
 *             will be null-terminated).
 * @param[in] bufsz The size of the buf string (i.e. number of chars). This 
 *             should be at least 4.
 * @return If a value other than zero is returned, then we couldn't match 
 *         an (int) id to the given planet.
 */
int planet_to_naif_idstr(dso::Planet p, char *buf,
                                    int bufsz) noexcept;

/* @brief Load a given cspice kernel if not already loaded.
 * @see
 * https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/kernel.html#Kernel%20Types
 */
int load_if_unloaded(const char *kernel) noexcept;

/* @brief Compute Ephemeris Time from TT passed as Julian Date via cspice
 * @note The function expects that:
 *       * An lsk (aka leap-second) kernel is already loaded, so that we
 *       can perform datetime computations
 */
inline double mjdtt2et(const dso::MjdEpoch &tt) noexcept {
  /* get ephemeris (ET) time from TT; assumes an LSK kernel is loaded ... */
  return unitim_c(tt.julian_date(), "JDTDT", "ET");
}

/* @brief Position of a target body relative to an observing body.
 * 
 * The reference frame for the returned vector is J2000 and the (returned)
 * values have units [km]. Note that no aberration corrections are applied.
 * @note The function expects that:
 *       * An spk kernel is already loaded so that we can compute the
 *         target/oberver positions
 * For the target/oberver ID's, see
 * https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html
 * Note that according to CSPICE documentation, J20i00 frame is actually ICRF,
 * see
 * https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/Tutorials/pdf/individual_docs/04_concepts.pdf
 *
 * @param[in] et The epoch as ephemeris time (see mjdtt2et)
 * @param[in] target_id An int representing the target body, using the
 *            interbal scpice IDs (see planet_to_naif_id)
 * @param[in] observer_id An int representing the observing body, using the
 *            interbal scpice IDs (planet_to_naif_id)
 * @param[out] pos A vector of size >= 3, where the X/Y/Z components of the
 *            vector are stored; units [km]
 */
inline int j2planet_pos_from(double et, int target_id, int observer_id,
                             double *pos) noexcept {
  double dummy;
  /* get the position of the planet */
  spkezp_c(target_id, et, "J2000", "NONE", observer_id, pos, &dummy);
  return 0;
}

/* @brief State of a target body relative to an observing body.
 * 
 * The reference frame for the returned vector is J2000 and the (returned)
 * values have units [km] and [km/sec]. Note that no aberration corrections 
 * are applied.
 *
 * @note The function expects that:
 *       * An spk kernel is already loaded so that we can compute the
 *         target/oberver positions
 * 
 * For the target/oberver ID's, see
 * https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html
 * Note that according to CSPICE documentation, J20i00 frame is actually ICRF,
 * see
 * https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/Tutorials/pdf/individual_docs/04_concepts.pdf
 *
 * @param[in] et The epoch as ephemeris time (see mjdtt2et)
 * @param[in] target_id An int as string representing the target body, using 
 *            the internal scpice IDs (see planet_to_naif_idstr). E.g. for 
 *            Earth, which has an id of 399, this should be "399".
 * @param[in] observer_id An int as string representing the observing body, 
 *            using the internal scpice IDs (planet_to_naif_id)
 * @param[out] pv A vector of size >= 6, where the X/Y/Z components of the
 *            position and velocity vectors are stored; units are [km] and 
 *            [km/sec].
 */
inline int j2planet_state_from(double et, const char *target_id,
                               const char *observer_id, double *pv) noexcept {
  double dummy;
  /* get the state of the planet */
  spkezr_c(target_id, et, "J2000", "NONE", observer_id, pv, &dummy);
  return 0;
}

} /* namespace cspice */

/* @brief Load a cspice kernel (of any type)
 * @param[in] kernel The name of the cspice kernel, see
 * https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/kernel.html#Kernel%20Types
 * @note This is a wrapper function around furnsh_c; see
 * https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/furnsh_c.html
 */
inline int load_spice_kernel(const char *kernel) noexcept {
  cspice::load_if_unloaded(kernel);
  return 0;
}

/* @brief Get planet position w.r.t Earth in GCRF/J2000
 * @param[in] p A Planet enum, denoting the target planet
 * @param[in] mjd_tt The date of request as an MJD in TT
 * @param[out] pos The poistion of the (target) planet at the date of request
 *            in [m] at the GCRF/J2000 frame
 */
int planet_pos(Planet p, const MjdEpoch &mjd_tt,
               Eigen::Matrix<double, 3, 1> &pos) noexcept;

/* @brief Get planet state (position+velocity) in GCRF/J2000
 * @param[in] p A Planet enum, denoting the target planet
 * @param[in] mjd_tt The date of request as an MJD in TT
 * @param[out] pv The state vector of the (target) planet at the date of 
 *            request in [m] and [m/sec] at the GCRF/J2000 frame.
 */
int planet_state(Planet p, const MjdEpoch &mjd_tt,
               Eigen::Matrix<double, 6, 1> &pv) noexcept;

/* @brief Get the Sun's and Moon's standard gravitational parameter (μ) off
 *        from a CSPICE/NAIF PCK Kernel
 * @param[in] pck_kernel A PCK kernel filename. If the kernel is already
 *         loaded, we will get the values without reloading it. If the
 *         filename is NULL, then we will try retrieving the constants
 *         assuming a PCK kernel is already loaded.
 * @param[out] GMSun Sun's gravitational constant according to the kernel in
 *         [km^3 / sec^2]
 * @param[out] GMMoon Moon's gravitational constant according to the kernel in
 *         [km^3 / sec^2]
 * @param[in] use_si Transform Sun's and Moon's gravitational constant to SI
 *         units, i.e. to m^3/sec^2 (instead of km^3/sec^2)
 * @return Anything other than zero denotes an error.
 */
int sun_moon_gm(double &GMSun, double &GMMoon, int use_si = false,
                const char *pck_kernel = nullptr) noexcept;

} /* namespace dso */

#endif
