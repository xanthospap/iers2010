#include "planets.hpp"
#include <algorithm>
#include <array>
#include <limits>

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
