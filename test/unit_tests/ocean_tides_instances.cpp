#include "ocean_tide.hpp"
#include <cstring>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

constexpr const int m1_deg = 20;
constexpr const int m1_ord = 16;
constexpr const int m2_deg = 17;
constexpr const int m2_ord = 12;
constexpr const int m3_deg = 12;
constexpr const int m3_ord = 12;

char *make_fn(const char *parta, const char *partb, char *buf) {
  std::strcpy(buf, parta);
  std::strcat(buf, "/");
  std::strcat(buf, partb);
  return buf;
}

int main([[maybe_unused]] int argc, char *argv[]) {
  char bf1[512], bf2[512], bf3[512], bf4[512];

  // Model1, no admittance
  {
    dso::OceanTide o1 = dso::groops_ocean_atlas(
        make_fn(argv[1], "model1/model1_001fileList.txt", bf1),
        make_fn(argv[1], "model1/", bf2));
    int max_order;
    assert(o1.max_degree(max_order) == m1_deg);
    assert(max_order == m1_ord);
  }
  
  // Model2, with admittance
  {
    dso::OceanTide o1 = dso::groops_ocean_atlas(
        make_fn(argv[1], "model2/model2_001fileList.txt", bf1),
        make_fn(argv[1], "model2/model2_002doodson.txt", bf2),
        make_fn(argv[1], "model2/model2_003admittance.txt", bf3),
        make_fn(argv[1], "model2/", bf4));
    int max_order;
    assert(o1.max_degree(max_order) == m2_deg);
    assert(max_order == m2_ord);
  }
  
  // Model3, with admittance
  {
    dso::OceanTide o1 = dso::groops_ocean_atlas(
        make_fn(argv[1], "model3/model3_001fileList.txt", bf1),
        make_fn(argv[1], "model3/model3_002doodson.txt", bf2),
        make_fn(argv[1], "model3/model3_003admittance.txt", bf3),
        make_fn(argv[1], "model3/", bf4));
    int max_order;
    assert(o1.max_degree(max_order) == m3_deg);
    assert(max_order == m3_ord);
  }

  return 0;
}
