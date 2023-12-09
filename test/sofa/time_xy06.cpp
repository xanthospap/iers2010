#include "geodesy/units.hpp"
#include "iau.hpp"
#include "sofa.h"
#include <cassert>
#include <chrono>
#include <cstdio>
#include <datetime/calendar.hpp>
#include <random>
#include <unistd.h>
#include <cstring>

using namespace std::chrono;

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<> id(0, 1000);
std::uniform_real_distribution<double> dd(-1e0, 1e0);

struct Clock {
  long mruntime;
  long msamples;
  time_point<high_resolution_clock> mstart;
  time_point<high_resolution_clock> mstop;
  char mname[24] = {'\0'};

  Clock(const char *n) { std::strcpy(mname, n); }
  Clock(const Clock &other) {
    mruntime = other.mruntime;
    msamples = other.msamples;
    mstart = other.mstart;
    mstop = other.mstop;
    std::strcpy(mname, other.mname);
  }
  Clock(Clock &&other) {
    mruntime = other.mruntime;
    msamples = other.msamples;
    mstart = other.mstart;
    mstop = other.mstop;
    std::strcpy(mname, other.mname);
  }
  Clock &operator=(const Clock &other) {
    if (this != &other) {
      mruntime = other.mruntime;
      msamples = other.msamples;
      mstart = other.mstart;
      mstop = other.mstop;
      std::strcpy(mname, other.mname);
    }
    return *this;
  }
  void start() { mstart = high_resolution_clock::now(); }
  auto stop() { mstop = high_resolution_clock::now(); }
  
  auto add_sample() {
    auto duration = duration_cast<microseconds>(mstop - mstart).count();
    mruntime += duration;
    mruntime /= 2;
  }
  
  bool is_faster(const Clock &other) const noexcept {
    return mruntime < other.mruntime;
  }
  
  void compare(const Clock &other) const noexcept {
    double ref = (double)(other.mruntime);
    double y = (double)mruntime;
    const double dx = 1e2 * (y-ref) / ref;
    printf("%s wrt %s : %+.1f %%\n", mname, other.mname, dx);
  }
};

constexpr const int num_tests = 50'00;
constexpr const double MAX_ARCSEC = 1e-9;
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
  c.add_sample();
  // printf("SOFA run-time = %15ld\n", duration.count());
  return dummy;
}

int run_mine1_timer(const std::vector<dso::MjdEpoch> &epochs, Clock &c) {
  assert(epochs.size() == num_tests);
  int dummy = 0;
  c.start();
  for (int i = 0; i < num_tests; i++) {
    /* check MINE */
    double xcip, ycip;
    dso::xycip06a(epochs[i], xcip, ycip);
    dummy -= (int)(xcip + ycip);
    /* flush cache */
    dummy += flush_cache(i < 200 ? i : 200);
  }
  c.stop();
  c.add_sample();
  //printf("MINE run-time = %15ld [v.normal]\n", duration.count());
  return dummy;
}

int run_mine2_timer(const std::vector<dso::MjdEpoch> &epochs, Clock &c) {
  assert(epochs.size() == num_tests);
  int dummy = 0;
  c.start();
  for (int i = 0; i < num_tests; i++) {
    /* check MINE */
    double xcip, ycip;
    dso::xycip06a_thread(epochs[i], xcip, ycip);
    dummy -= (int)(xcip + ycip);
    /* flush cache */
    dummy += flush_cache(i < 200 ? i : 200);
  }
  c.stop();
  c.add_sample();
  //printf("MINE run-time = %15ld [v.thread]\n", duration.count());
  return dummy;
}

int run_mine3_timer(const std::vector<dso::MjdEpoch> &epochs, Clock &c) {
  assert(epochs.size() == num_tests);
  int dummy = 0;
  c.start();
  for (int i = 0; i < num_tests; i++) {
    /* check MINE */
    double xcip, ycip;
    dso::xycip06a_sofa(epochs[i], xcip, ycip);
    dummy -= (int)(xcip + ycip);
    /* flush cache */
    dummy += flush_cache(i < 200 ? i : 200);
  }
  c.stop();
  c.add_sample();
  //printf("MINE run-time = %15ld [v.sofa]\n", duration.count());
  return dummy;
}

int run_mine4_timer(const std::vector<dso::MjdEpoch> &epochs, Clock &c) {
  assert(epochs.size() == num_tests);
  int dummy = 0;
  c.start();
  for (int i = 0; i < num_tests; i++) {
    /* check MINE */
    double xcip, ycip;
    dso::xycip06a_new(epochs[i], xcip, ycip);
    dummy -= (int)(xcip + ycip);
    /* flush cache */
    dummy += flush_cache(i < 200 ? i : 200);
  }
  c.stop();
  c.add_sample();
  //printf("MINE run-time = %15ld [v.new]\n", duration.count());
  return dummy;
}

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

  int dummy;
  Clock c1("v1"), c2("v1.threads"), c3("sofa-copy"), c4("fusion"),
      c5("sofa-c");

  dummy += run_mine4_timer(mjd_epochs,c4);
  sleep(1);
  dummy += run_mine1_timer(mjd_epochs,c1);
  sleep(1);
  dummy += run_sofa_timer(jd_epochs,c5);
  sleep(1);
  dummy += run_mine2_timer(mjd_epochs,c2);
  sleep(1);
  dummy += run_mine3_timer(mjd_epochs,c3);
  sleep(1);

  dummy += run_sofa_timer(jd_epochs,c5);
  sleep(1);
  dummy += run_mine4_timer(mjd_epochs,c4);
  sleep(1);
  dummy += run_mine2_timer(mjd_epochs,c2);
  sleep(1);
  dummy += run_mine3_timer(mjd_epochs,c3);
  sleep(1);
  dummy += run_mine1_timer(mjd_epochs,c1);
  sleep(1);

  dummy += run_mine2_timer(mjd_epochs,c2);
  sleep(1);
  dummy += run_mine3_timer(mjd_epochs,c3);
  sleep(1);
  dummy += run_mine4_timer(mjd_epochs,c4);
  sleep(1);
  dummy += run_mine1_timer(mjd_epochs,c1);
  sleep(1);
  dummy += run_sofa_timer(jd_epochs,c5);
  sleep(1);

  dummy += run_mine3_timer(mjd_epochs,c3);
  sleep(1);
  dummy += run_mine1_timer(mjd_epochs,c1);
  sleep(1);
  dummy += run_sofa_timer(jd_epochs,c5);
  sleep(1);
  dummy += run_mine4_timer(mjd_epochs,c4);
  sleep(1);
  dummy += run_mine2_timer(mjd_epochs,c2);
  sleep(1);

  dummy += run_mine1_timer(mjd_epochs,c1);
  sleep(1);
  dummy += run_sofa_timer(jd_epochs,c5);
  sleep(1);
  dummy += run_mine2_timer(mjd_epochs,c2);
  sleep(1);
  dummy += run_mine3_timer(mjd_epochs,c3);
  sleep(1);
  dummy += run_mine4_timer(mjd_epochs,c4);
  sleep(1);

  // choose faster clock
  Clock ref(c1);
  if (! ref.is_faster(c2) ) ref = c2;
  if (! ref.is_faster(c3) ) ref = c3;
  if (! ref.is_faster(c4) ) ref = c4;
  if (! ref.is_faster(c5) ) ref = c5;
  c1.compare(ref);
  c2.compare(ref);
  c4.compare(ref);
  c5.compare(ref);

  return 0;
}
