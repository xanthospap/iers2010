#include "ocean_tide.hpp"
#include <cstring>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <random>
#include <cmath>

constexpr const int m1_deg = 20;
constexpr const int m1_ord = 16;
constexpr const int m2_deg = 17;
constexpr const int m2_ord = 12;
constexpr const int m3_deg = 12;
constexpr const int m3_ord = 12;

std::uniform_real_distribution<double> unif(0e0, 2*M_PI);
std::default_random_engine re;

char *make_fn(const char *parta, const char *partb, char *buf) {
  std::strcpy(buf, parta);
  std::strcat(buf, "/");
  std::strcat(buf, partb);
  return buf;
}

int main([[maybe_unused]] int argc, char *argv[]) {
  char bf1[512], bf2[512], bf3[512], bf4[512];
  dso::MjdEpoch t = dso::MjdEpoch::random(dso::modified_julian_day(53005),
                                          dso::modified_julian_day(60310));
  double delaunay[6];
  for (int i=0; i<6; i++) delaunay[i] = unif(re);

  // Model1, no admittance
  {
    dso::OceanTide o1 = dso::groops_ocean_atlas(
        make_fn(argv[1], "model1/model1_001fileList.txt", bf1),
        make_fn(argv[1], "model1/", bf2));
    int max_order;
    assert(o1.atlas().max_atlas_degree(max_order) == m1_deg);
    assert(max_order == m1_ord);
    assert(o1.max_degree() == m1_deg);
    assert(o1.max_order() == m1_ord);
  }
  
  // Model2, with admittance
  {
    dso::OceanTide o1 = dso::groops_ocean_atlas(
        make_fn(argv[1], "model2/model2_001fileList.txt", bf1),
        make_fn(argv[1], "model2/model2_002doodson.txt", bf2),
        make_fn(argv[1], "model2/model2_003admittance.txt", bf3),
        make_fn(argv[1], "model2/", bf4));
    int max_order;
    assert(o1.atlas().max_atlas_degree(max_order) == m2_deg);
    assert(max_order == m2_ord);
    assert(o1.max_degree() == m2_deg);
    assert(o1.max_order() == m2_ord);
  }
  
  // Model3, with admittance
  {
    dso::OceanTide o1 = dso::groops_ocean_atlas(
        make_fn(argv[1], "model3/model3_001fileList.txt", bf1),
        make_fn(argv[1], "model3/model3_002doodson.txt", bf2),
        make_fn(argv[1], "model3/model3_003admittance.txt", bf3),
        make_fn(argv[1], "model3/", bf4));
    int max_order;
    assert(o1.atlas().max_atlas_degree(max_order) == m3_deg);
    assert(max_order == m3_ord);
    assert(o1.max_degree() == m3_deg);
    assert(o1.max_order() == m3_ord);
  }

  // Copy Con'tors
  {
    dso::OceanTide o1 = dso::groops_ocean_atlas(
        make_fn(argv[1], "model3/model3_001fileList.txt", bf1),
        make_fn(argv[1], "model3/model3_002doodson.txt", bf2),
        make_fn(argv[1], "model3/model3_003admittance.txt", bf3),
        make_fn(argv[1], "model3/", bf4));
    int max_order;
    assert(o1.atlas().max_atlas_degree(max_order) == m3_deg);
    assert(max_order == m3_ord);
    assert(o1.max_degree() == m3_deg);
    assert(o1.max_order() == m3_ord);

    dso::OceanTide o2(o1);
    int max_order2;
    assert(o2.atlas().max_atlas_degree(max_order2) == m3_deg);
    assert(max_order2 == m3_ord);
    assert(o2.max_degree() == m3_deg);
    assert(o2.max_order() == m3_ord);

    const double dut1 = unif(re);
    assert(!o1.stokes_coeffs(t, t.tt2ut1(dut1), delaunay));
    assert(!o2.stokes_coeffs(t, t.tt2ut1(dut1), delaunay));

    assert(o1.stokes_coeffs().max_degree() == o2.stokes_coeffs().max_degree());
    assert(o1.stokes_coeffs().max_order() == o2.stokes_coeffs().max_order());
    for (int i = 0; i <= o1.stokes_coeffs().max_degree(); i++) {
      for (int j = 0; j <= std::min(i, o1.stokes_coeffs().max_order()); j++) {
        assert(std::abs(o1.stokes_coeffs().C(i, j) -
                        o2.stokes_coeffs().C(i, j)) < 1e-15);
        assert(std::abs(o1.stokes_coeffs().S(i, j) -
                        o2.stokes_coeffs().S(i, j)) < 1e-15);
      }
    }
  }
  
  // Assignment
  {
    dso::OceanTide o3 = dso::groops_ocean_atlas(
        make_fn(argv[1], "model3/model3_001fileList.txt", bf1),
        make_fn(argv[1], "model3/model3_002doodson.txt", bf2),
        make_fn(argv[1], "model3/model3_003admittance.txt", bf3),
        make_fn(argv[1], "model3/", bf4));
    int max_order;
    assert(o3.atlas().max_atlas_degree(max_order) == m3_deg);
    assert(max_order == m3_ord);
    assert(o3.max_degree() == m3_deg);
    assert(o3.max_order() == m3_ord);

    dso::OceanTide o2 = o3;
    int max_order2;
    assert(o2.atlas().max_atlas_degree(max_order2) == m3_deg);
    assert(max_order2 == m3_ord);
    assert(o2.max_degree() == m3_deg);
    assert(o2.max_order() == m3_ord);
    
    dso::OceanTide o1 = dso::groops_ocean_atlas(
        make_fn(argv[1], "model1/model1_001fileList.txt", bf1),
        make_fn(argv[1], "model1/", bf2));
    assert(o1.atlas().max_atlas_degree(max_order) == m1_deg);
    assert(max_order == m1_ord);
    assert(o1.max_degree() == m1_deg);
    assert(o1.max_order() == m1_ord);
    
    o2 = o1;
    assert(o2.atlas().max_atlas_degree(max_order2) == m1_deg);
    assert(max_order2 == m1_ord);
    assert(o2.max_degree() == m1_deg);
    assert(o2.max_order() == m1_ord);

    const double dut1 = unif(re);
    assert(!o1.stokes_coeffs(t, t.tt2ut1(dut1), delaunay));
    assert(!o2.stokes_coeffs(t, t.tt2ut1(dut1), delaunay));

    assert(o1.stokes_coeffs().max_degree() == o2.stokes_coeffs().max_degree());
    assert(o1.stokes_coeffs().max_order() == o2.stokes_coeffs().max_order());
    for (int i = 0; i <= o1.stokes_coeffs().max_degree(); i++) {
      for (int j = 0; j <= std::min(i, o1.stokes_coeffs().max_order()); j++) {
        assert(std::abs(o1.stokes_coeffs().C(i, j) -
                        o2.stokes_coeffs().C(i, j)) < 1e-15);
        assert(std::abs(o1.stokes_coeffs().S(i, j) -
                        o2.stokes_coeffs().S(i, j)) < 1e-15);
      }
    }
  }

  return 0;
}
