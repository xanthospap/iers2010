#include "hardisp.hpp"
#include <cassert>
#include <cstring>
#include <algorithm>

// extracted from test file NTUA.BLQ
const dso::BlqSiteInfo bdyng{
    "DYNG",
    {.00403, .00159, .00094, .00042, .00065, .00043, .00021, .00011, .00023,
     .00010, .00008, .00136, .00032, .00029, .00009, .00039, .00020, .00013,
     .00002, .00003, .00002, .00002, .00074, .00031, .00016, .00009, .00048,
     .00020, .00016, .00004, .00004, .00003, .00002},
    {dso::deg2rad(-86.8),  dso::deg2rad(-53.5),  dso::deg2rad(-112.3),
     dso::deg2rad(-65.3),  dso::deg2rad(-98.8),  dso::deg2rad(-149.2),
     dso::deg2rad(-102.7), dso::deg2rad(-135.0), dso::deg2rad(-12.1),
     dso::deg2rad(-2.8),   dso::deg2rad(-1.0),   dso::deg2rad(55.0),
     dso::deg2rad(63.1),   dso::deg2rad(41.7),   dso::deg2rad(57.9),
     dso::deg2rad(130.3),  dso::deg2rad(90.0),   dso::deg2rad(128.1),
     dso::deg2rad(110.3),  dso::deg2rad(-163.5), dso::deg2rad(-159.8),
     dso::deg2rad(-176.1), dso::deg2rad(-89.4),  dso::deg2rad(-53.9),
     dso::deg2rad(-121.5), dso::deg2rad(-72.2),  dso::deg2rad(38.0),
     dso::deg2rad(44.6),   dso::deg2rad(38.0),   dso::deg2rad(55.4),
     dso::deg2rad(-141.0), dso::deg2rad(-162.5), dso::deg2rad(-174.5)}};
const dso::BlqSiteInfo bzywi{
    "ZYWI",
    {.00405, .00143, .00088, .00036, .00168, .00096, .00056, .00011, .00049,
     .00027, .00022, .00149, .00033, .00035, .00009, .00039, .00030, .00013,
     .00005, .00005, .00003, .00003, .00033, .00012, .00005, .00003, .00035,
     .00012, .00012, .00001, .00002, .00001, .00000},
    {dso::deg2rad(-72.0),  dso::deg2rad(-39.6),  dso::deg2rad(-92.0),
     dso::deg2rad(-50.3),  dso::deg2rad(-62.3),  dso::deg2rad(-101.4),
     dso::deg2rad(-63.1),  dso::deg2rad(-144.1), dso::deg2rad(8.8),
     dso::deg2rad(5.8),    dso::deg2rad(1.1),    dso::deg2rad(69.4),
     dso::deg2rad(104.8),  dso::deg2rad(44.1),   dso::deg2rad(94.6),
     dso::deg2rad(105.9),  dso::deg2rad(43.0),   dso::deg2rad(104.3),
     dso::deg2rad(3.1),    dso::deg2rad(-169.8), dso::deg2rad(-166.2),
     dso::deg2rad(-177.2), dso::deg2rad(-78.1),  dso::deg2rad(-35.3),
     dso::deg2rad(-109.8), dso::deg2rad(-55.6),  dso::deg2rad(39.9),
     dso::deg2rad(4.9),    dso::deg2rad(39.5),   dso::deg2rad(63.9),
     dso::deg2rad(-171.0), dso::deg2rad(168.1),  dso::deg2rad(-174.9)}};

bool site_blq_equal(const dso::BlqSiteInfo &b1, const dso::BlqSiteInfo &b2) {
  if (std::strcmp(b1.site, b2.site)) return false;
  for (int i=0; i<11*3; i++) {
    if (b1.amplitudes[i] != b2.amplitudes[i])
      return false;
  }
  for (int i=0; i<11*3; i++) {
    if (b1.phases[i] != b2.phases[i])
      return false;
  }
  return true;
}

int main(int argc, char *argv[]) {
  if (argc!=2) {
    fprintf(stderr, "Usage: %s [BLQ FILE]\n", argv[0]);
    return 1;
  }

  printf("Testing reading of BLQ file\n");

  std::vector<dso::BlqSiteInfo> blqInfoVec;

  // read all entries off from the BLQ file
  if (dso::parse_blq(argv[1], blqInfoVec, nullptr)) {
    fprintf(stderr, "Failed parsing BLQ file\n");
    return 1;
  }

  // check number of sites collected
  assert(blqInfoVec.size() == 357);

  // we should have an entry for DYNG
  auto it = std::find_if(blqInfoVec.begin(), blqInfoVec.end(),
                         [](const dso::BlqSiteInfo &b) {
                           return !std::strcmp(b.site, "DYNG");
                         });
  assert( it != blqInfoVec.end());
  assert( site_blq_equal(*it, bdyng) );
  
  // we should have an entry for ZYWI
 it = std::find_if(blqInfoVec.begin(), blqInfoVec.end(),
                         [](const dso::BlqSiteInfo &b) {
                           return !std::strcmp(b.site, "ZYWI");
                         });
  assert( it != blqInfoVec.end());
  assert( site_blq_equal(*it, bzywi) );

  // let's extract only two sites
  const char *dyng = "DYNG";
  const char *zywi = "ZYWI";
  std::vector<const char *> sites;
  sites.emplace_back(dyng);
  sites.emplace_back(zywi);
  // read all entries off from the BLQ file
  if (dso::parse_blq(argv[1], blqInfoVec, &sites)) {
    fprintf(stderr, "Failed parsing BLQ file\n");
    return 1;
  }

  // check number of sites collected
  assert(blqInfoVec.size() == 2);

  // we should have an entry for DYNG
 it = std::find_if(blqInfoVec.begin(), blqInfoVec.end(),
                         [](const dso::BlqSiteInfo &b) {
                           return !std::strcmp(b.site, "DYNG");
                         });
  assert( it != blqInfoVec.end());
  assert( site_blq_equal(*it, bdyng) );
  
  // we should have an entry for ZYWI
 it = std::find_if(blqInfoVec.begin(), blqInfoVec.end(),
                         [](const dso::BlqSiteInfo &b) {
                           return !std::strcmp(b.site, "ZYWI");
                         });
  assert( it != blqInfoVec.end());
  assert( site_blq_equal(*it, bzywi) );

  // if we pass in an empty vector, we should get an empty vector back
  sites.clear();
  if (dso::parse_blq(argv[1], blqInfoVec, &sites)) {
    fprintf(stderr, "Failed parsing BLQ file\n");
    return 1;
  }
  assert( !blqInfoVec.size() );
  
  // let's try for a site that isn't there
  const char *miss = "NONE";
  sites.emplace_back(dyng);
  sites.emplace_back(miss);
  if (dso::parse_blq(argv[1], blqInfoVec, &sites)) {
    fprintf(stderr, "Failed parsing BLQ file\n");
    return 1;
  }
  assert( blqInfoVec.size() == 1);
  assert( site_blq_equal(*(blqInfoVec.begin()), bdyng) );

  printf("Test ok\n");

  return 0;
}
