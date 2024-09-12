#include "tide_atlas.hpp"
#include <cassert>
#include <cstdio>
#include <algorithm>
#include <cstring>

using namespace dso;

int main() {

  char buf[64]={"\0"};
  const char *wstr[] = {"011.111", "111.111", "211.111", "311.111", "411.111"};
  const char *sstr[] = {"s1", "s2", "Om1", "foo", "bar1"};
  const int dszs[] = {20,30,40,50,60};
  const int oszs[] = {20,20,30,40,50};

  TideAtlas atlas;
  atlas.waves().reserve(5);
  atlas.append_wave(TidalWave(DoodsonConstituent(0,1,1,1,1,1,0),2,2,"s1"  ),20,20);
  atlas.append_wave(TidalWave(DoodsonConstituent(1,1,1,1,1,1,0),2,2,"s2"  ),30,20);
  atlas.append_wave(TidalWave(DoodsonConstituent(2,1,1,1,1,1,0),2,2,"Om1" ),40,30);
  atlas.append_wave(TidalWave(DoodsonConstituent(3,1,1,1,1,1,0),2,2,"foo" ),50,40);
  atlas.append_wave(TidalWave(DoodsonConstituent(4,1,1,1,1,1,0),2,2,"bar1"),60,50);

  /* this should throw ! */
  try {
    atlas.append_wave(TidalWave(DoodsonConstituent(0,1,1,1,1,1,0)),20,20);
    assert(2==1);
  } catch (std::exception &) {
    fprintf(stderr, "Exception thrown and caught; its all good, expected!\n");
    ;
  }

  /* assert size of atlas */
  assert(atlas.waves().size() == 5);

  int j=0;
  for (auto it = atlas.waves().begin(); it != atlas.waves().end(); ++it) {
    //printf("Wave <%s> aka <", it->wave().doodson().str(buf, false));
    //for (int i=0 ;i<6; i++) printf("%d,", it->wave().doodson().int_array()[i]);
    //printf(">\n");
    assert(!std::strcmp(it->wave().doodson().str(buf, false), wstr[j]));
    assert(!std::strcmp(it->name(), sstr[j]));
    assert((it->stokes_sin().max_degree() == it->stokes_cos().max_degree()) &&
           (it->stokes_cos().max_degree() == dszs[j]));
    assert((it->stokes_sin().max_order() == it->stokes_cos().max_order()) &&
           (it->stokes_cos().max_order() == oszs[j]));
    ++j;
  }

  auto it = atlas.find_tidal_wave(DoodsonConstituent(3,1,1,1,1,1,0));
  assert(it == atlas.waves().begin() + 3);
  
  atlas.append_wave(TidalWave(DoodsonConstituent(3,1,1,1,1,2,0),2,2,"foo" ),50,40);
  atlas.append_wave(TidalWave(DoodsonConstituent(4,1,1,1,1,2,0),2,2,"bar1"),60,50);
  assert(atlas.waves().size() == 7);
  it = atlas.find_tidal_wave(DoodsonConstituent(3,1,1,1,1,1,0));
  assert(it == atlas.waves().begin() + 3);

  return 0;
}
