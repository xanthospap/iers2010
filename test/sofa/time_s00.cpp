#include "iau.hpp"
#include "sofa.h"
#include "geodesy/units.hpp"
#include <cstdio>
#include <cassert>
#include <datetime/calendar.hpp>
#include <random>
#include <chrono>
using namespace std::chrono;

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<> id(0, 1000);
std::uniform_real_distribution<double > dd(-1e0, 1e0);

constexpr const int num_tests = 10'00;
constexpr const double MAX_ARCSEC = 1e-9;
struct JD {
  double jd1;
  double jd2;
  JD(double d1, double d2): jd1(d1), jd2(d2) {};
};

int flush_cache(int n) noexcept {
  if (n<200) n = 200;
  double *ar = new double[n];
  for (int i=0; i<n; i++) ar[i] = dd(gen);
  int d = 0;
  for (int i=0; i<n; i++) d += ar[i];
  delete[] ar;
  return d >= 0 ? 1 : 0;
}

int run_sofa_timer(const std::vector<JD> &epochs) {
  assert(epochs.size() == num_tests);
  int dummy = 0;
  auto start = high_resolution_clock::now();
  for (int i = 0; i < num_tests; i++) {
    /* Two-part JD */
    const double dj1 = epochs[i].jd1;
    const double dj2 = epochs[i].jd2;
    /* check SOFA */
    const auto s = iauS00(dj1, dj2, 0, 0);
    dummy -= (int)s;
    /* flush cache */
    dummy += flush_cache(i < 200 ? i : 200);
  }
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  printf("SOFA run-time = %15ld\n", duration.count());
  return dummy;
}

int run_mine_timer(const std::vector<dso::MjdEpoch> &epochs) {
  assert(epochs.size() == num_tests);
  int dummy = 0;
  auto start = high_resolution_clock::now();
  for (int i = 0; i < num_tests; i++) {
    /* check MINE */
    const auto s = dso::s06(epochs[i], 0, 0);
    dummy -= (int)s;
    /* flush cache */
    dummy += flush_cache(i < 200 ? i : 200);
  }
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  printf("MINE run-time = %15ld\n", duration.count());
  return dummy;
}

//int run_thrd_timer(const std::vector<dso::MjdEpoch> &epochs) {
//  assert(epochs.size() == num_tests);
//  int dummy = 0;
//  auto start = high_resolution_clock::now();
//  for (int i = 0; i < num_tests; i++) {
//    /* check MINE */
//    double xcip, ycip, s;
//    dso::xys06a(epochs[i], xcip, ycip, s);
//    dummy -= (int)(s-xcip+ycip);
//    /* flush cache */
//    dummy += flush_cache(i < 200 ? i : 200);
//  }
//  auto stop = high_resolution_clock::now();
//  auto duration = duration_cast<microseconds>(stop - start);
//  printf("MINE run-time = %15ld\n", duration.count());
//  return dummy;
//}

int main() {

  /* generate random but corresponing dates */
  std::vector<dso::MjdEpoch> mjd_epochs;
  for (int i = 0; i < num_tests; i++) {
    mjd_epochs.emplace_back(dso::MjdEpoch::random(40587, 66154));
  }
  std::vector<JD> jd_epochs;
  for (int i = 0; i < num_tests; i++) {
    jd_epochs.emplace_back(mjd_epochs[i].imjd() + dso::MJD0_JD,
                            mjd_epochs[i].fractional_days());
  }

    /* time SOFA */
    run_sofa_timer(jd_epochs);

    /* time this implementation */
    run_mine_timer(mjd_epochs);
    
    /* time multi-thread implementation */
    //run_thrd_timer(mjd_epochs);
    
    /* time SOFA */
    run_sofa_timer(jd_epochs);

    /* time this implementation */
    run_mine_timer(mjd_epochs);
    
    /* time multi-thread implementation */
    //run_thrd_timer(mjd_epochs);
    
    /* time SOFA */
    run_sofa_timer(jd_epochs);

    /* time this implementation */
    run_mine_timer(mjd_epochs);
    
    /* time multi-thread implementation */
    //run_thrd_timer(mjd_epochs);

  return 0;
}
