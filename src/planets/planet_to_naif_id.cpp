#include "planets.hpp"
#include <algorithm>
#include <array>
#include <limits>
#include <charconv>
#include <cstring>

namespace {
/* Match Planet enums to CSPICE IDs, see
 * https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html
 */
struct Planet2NaifId {
  dso::Planet p;
  int spiceId;
};
constexpr const std::array<Planet2NaifId, 8> PlanetId = {
    {{dso::Planet::EARTH, 399},
     {dso::Planet::MOON, 301},
     {dso::Planet::SUN, 10},
     {dso::Planet::MERCURY, 199},
     {dso::Planet::VENUS, 299},
     {dso::Planet::MARS, 499},
     {dso::Planet::JUPITER, 599},
     {dso::Planet::SATURN, 699}}};
} /* unnamed namespace */

int dso::cspice::planet_to_naif_id(dso::Planet p, int &id) noexcept {
  auto it = std::find_if(PlanetId.begin(), PlanetId.end(),
                         [&](const Planet2NaifId &pid) { return pid.p == p; });
  id = (it != PlanetId.end()) ? it->spiceId : std::numeric_limits<int>::max();

  return (id == std::numeric_limits<int>::max()) ? 1 : 0;
}

int dso::cspice::planet_to_naif_idstr(dso::Planet p, char *buf,
                                    int bufsz) noexcept {
  /* string to null */
  std::memset(buf, '\0', bufsz);
  /* get the planet's id */
  int id;
  if (dso::cspice::planet_to_naif_id(p, id))
    return 1;
  /* transform the integer value to a string and write it out to buf */
  if (auto [ptr, ec] = std::to_chars(buf, buf + bufsz, id); ec != std::errc())
    return 2;

  return 0;
}
