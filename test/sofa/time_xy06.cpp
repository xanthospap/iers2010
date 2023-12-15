#include "geodesy/units.hpp"
#include "iau.hpp"
#include "sofa.h"
#include <cassert>
#include <chrono>
#include <cstdio>
#include <random>
#include <unistd.h>
#include <cstring>
#include "clock.hpp"

using namespace std::chrono;

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<> id(0, 1000);
std::uniform_real_distribution<double> dd(-1e0, 1e0);

constexpr const int num_tests = 50'00;
// constexpr const double MAX_ARCSEC = 1e-9;
struct JD {
  double jd1;
  double jd2;
  JD(double d1, double d2) : jd1(d1), jd2(d2){};
};

int flush_cache(int n) noexcept {
  if (n < 200)
    n = 200;
  double *ar = new double[n];
  for (int i = 0; i < n; i++)
    ar[i] = dd(gen);
  int d = 0;
  for (int i = 0; i < n; i++)
    d += ar[i];
  delete[] ar;
  return d >= 0 ? 1 : 0;
}

int run_sofa_timer(const std::vector<JD> &epochs, Clock &c) {
  assert(epochs.size() == num_tests);
  int dummy = 0;
  c.start();
  for (int i = 0; i < num_tests; i++) {
    /* Two-part JD */
    const double dj1 = epochs[i].jd1;
    const double dj2 = epochs[i].jd2;
    /* check SOFA */
    double xcip, ycip;
    iauXy06(dj1, dj2, &xcip, &ycip);
    dummy -= (int)(xcip + ycip);
    /* flush cache */
    dummy += flush_cache(i < 200 ? i : 200);
  }
  c.stop();
  return dummy;
}

int run_mine1_timer(const std::vector<dso::MjdEpoch> &epochs, Clock &c) {
  assert(epochs.size() == num_tests);
  int dummy = 0;
  c.start();
  for (int i = 0; i < num_tests; i++) {
    double xcip, ycip;
    dso::xycip06a(epochs[i], xcip, ycip);
    dummy -= (int)(xcip + ycip);
    dummy += flush_cache(i < 200 ? i : 200);
  }
  c.stop();
  return dummy;
}

int run_mine2_timer(const std::vector<dso::MjdEpoch> &epochs, Clock &c) {
  assert(epochs.size() == num_tests);
  int dummy = 0;
  c.start();
  for (int i = 0; i < num_tests; i++) {
    double xcip, ycip;
    dso::extra::xycip06a(epochs[i], xcip, ycip);
    dummy -= (int)(xcip + ycip);
    dummy += flush_cache(i < 200 ? i : 200);
  }
  c.stop();
  return dummy;
}

//int run_mine3_timer(const std::vector<dso::MjdEpoch> &epochs, Clock &c) {
//  assert(epochs.size() == num_tests);
//  int dummy = 0;
//  c.start();
//  for (int i = 0; i < num_tests; i++) {
//    double xcip, ycip;
//    dso::xycip06a_sofa(epochs[i], xcip, ycip);
//    dummy -= (int)(xcip + ycip);
//    dummy += flush_cache(i < 200 ? i : 200);
//  }
//  c.stop();
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

  int dummy = 0;
  Clock c1("v1"), c2("v2"), c5("sofa-c"), cd("dummy");
  
  dummy += run_mine2_timer(mjd_epochs,c5); /* first call, no timer */
  sleep(1);

  dummy += run_mine1_timer(mjd_epochs,c1);
  sleep(1);
  dummy += run_sofa_timer(jd_epochs,c5);
  sleep(1);
  dummy += run_mine2_timer(mjd_epochs,c2);

  dummy += run_sofa_timer(jd_epochs,c5);
  sleep(1);
  dummy += run_mine2_timer(mjd_epochs,c2);
  sleep(1);
  dummy += run_mine1_timer(mjd_epochs,c1);
  sleep(1);

  dummy += run_mine2_timer(mjd_epochs,c2);
  sleep(1);
  dummy += run_mine1_timer(mjd_epochs,c1);
  sleep(1);
  dummy += run_sofa_timer(jd_epochs,c5);
  sleep(1);

  dummy += run_mine1_timer(mjd_epochs,c1);
  sleep(1);
  dummy += run_sofa_timer(jd_epochs,c5);
  sleep(1);
  dummy += run_mine2_timer(mjd_epochs,c2);
  sleep(1);

  dummy += run_mine1_timer(mjd_epochs,c1);
  sleep(1);
  dummy += run_sofa_timer(jd_epochs,c5);
  sleep(1);
  dummy += run_mine2_timer(mjd_epochs,c2);
  sleep(1);

  // choose faster clock
  Clock ref(c1);
  if (! ref.is_faster(c2) ) ref = c2;
  if (! ref.is_faster(c5) ) ref = c5;
  c1.compare(ref);
  c2.compare(ref);
  c5.compare(ref);

  // write dummy so that it is not optimized-out
  printf("%d\n", dummy);

  return 0;
}
