#include "hardisp.hpp"
#include <cassert>
#include <cstring>
#include <algorithm>

// extracted from test file NTUA.BLQ
const iers2010::BlqSiteInfo bdyng{
  "DYNG", 
  {.00403,.00159,.00094,.00042,.00065,.00043,.00021,.00011,.00023,.00010,.00008,
  .00136,.00032,.00029,.00009,.00039,.00020,.00013,.00002,.00003,.00002,.00002,
  .00074,.00031,.00016,.00009,.00048,.00020,.00016,.00004,.00004,.00003,.00002},
  {-86.8, -53.5,-112.3, -65.3, -98.8,-149.2,-102.7,-135.0, -12.1,  -2.8,  -1.0,
    55.0,  63.1,  41.7,  57.9, 130.3,  90.0, 128.1, 110.3,-163.5,-159.8,-176.1,
   -89.4, -53.9,-121.5, -72.2,  38.0,  44.6,  38.0,  55.4,-141.0,-162.5,-174.5}
};
const iers2010::BlqSiteInfo bzywi{
  "ZYWI",
  {.00405,.00143,.00088,.00036,.00168,.00096,.00056,.00011,.00049,.00027,.00022,
  .00149,.00033,.00035,.00009,.00039,.00030,.00013,.00005,.00005,.00003,.00003,
  .00033,.00012,.00005,.00003,.00035,.00012,.00012,.00001,.00002,.00001,.00000},
   {-72.0, -39.6, -92.0, -50.3, -62.3,-101.4, -63.1,-144.1,   8.8,   5.8,   1.1,
    69.4, 104.8,  44.1,  94.6, 105.9,  43.0, 104.3,   3.1,-169.8,-166.2,-177.2,
   -78.1, -35.3,-109.8, -55.6,  39.9,   4.9,  39.5,  63.9,-171.0, 168.1,-174.9}
};

bool site_blq_equal(const iers2010::BlqSiteInfo &b1, const iers2010::BlqSiteInfo &b2) {
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

  std::vector<iers2010::BlqSiteInfo> blqInfoVec;

  // read all entries off from the BLQ file
  if (iers2010::parse_blq(argv[1], blqInfoVec, nullptr)) {
    fprintf(stderr, "Failed parsing BLQ file\n");
    return 1;
  }

  // check number of sites collected
  assert(blqInfoVec.size() == 357);

  // we should have an entry for DYNG
  auto it = std::find_if(blqInfoVec.begin(), blqInfoVec.end(),
                         [](const iers2010::BlqSiteInfo &b) {
                           return !std::strcmp(b.site, "DYNG");
                         });
  assert( it != blqInfoVec.end());
  assert( site_blq_equal(*it, bdyng) );
  
  // we should have an entry for ZYWI
 it = std::find_if(blqInfoVec.begin(), blqInfoVec.end(),
                         [](const iers2010::BlqSiteInfo &b) {
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
  if (iers2010::parse_blq(argv[1], blqInfoVec, &sites)) {
    fprintf(stderr, "Failed parsing BLQ file\n");
    return 1;
  }

  // check number of sites collected
  assert(blqInfoVec.size() == 2);

  // we should have an entry for DYNG
 it = std::find_if(blqInfoVec.begin(), blqInfoVec.end(),
                         [](const iers2010::BlqSiteInfo &b) {
                           return !std::strcmp(b.site, "DYNG");
                         });
  assert( it != blqInfoVec.end());
  assert( site_blq_equal(*it, bdyng) );
  
  // we should have an entry for ZYWI
 it = std::find_if(blqInfoVec.begin(), blqInfoVec.end(),
                         [](const iers2010::BlqSiteInfo &b) {
                           return !std::strcmp(b.site, "ZYWI");
                         });
  assert( it != blqInfoVec.end());
  assert( site_blq_equal(*it, bzywi) );

  // if we pass in an empty vector, we should get an empty vector back
  sites.clear();
  if (iers2010::parse_blq(argv[1], blqInfoVec, &sites)) {
    fprintf(stderr, "Failed parsing BLQ file\n");
    return 1;
  }
  assert( !blqInfoVec.size() );
  
  // let's try for a site that isn't there
  const char *miss = "NONE";
  sites.emplace_back(dyng);
  sites.emplace_back(miss);
  if (iers2010::parse_blq(argv[1], blqInfoVec, &sites)) {
    fprintf(stderr, "Failed parsing BLQ file\n");
    return 1;
  }
  assert( blqInfoVec.size() == 1);
  assert( site_blq_equal(*(blqInfoVec.begin()), bdyng) );

  printf("Test ok\n");

  return 0;
}
