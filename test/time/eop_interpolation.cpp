#include "eop.hpp"
#include <cassert>
#include <cstdio>
#include "clock.hpp"
#include <unistd.h>
#include <random>

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<> id(0, 1000);
std::uniform_real_distribution<double> dd(-1e0, 1e0);

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

template<typename T>
int runimpl1(const T &eop, Clock &c) {
  c.start();
  int dummy = 0;
  double max;
  dso::MjdEpoch t(57752);
  while (t < dso::MjdEpoch(57754)) {
    dso::EopRecord ie;
    eop.interpolate(t, ie);
    max = (ie.xp() > ie.yp()) ? ie.xp() : ie.yp();
    max = (max > ie.dut()) ? max : ie.dut();
    max = (max > ie.lod()) ? max : ie.lod();
    max = (max > ie.dX()) ? max : ie.dX();
    t.add_seconds(30e0);
    /* flush cache */
    int i = (int)(max);
    dummy += flush_cache(i < 200 ? i : 200);
  }
  c.stop();
  c.add_sample();
  return dummy;
}

template<typename T>
int runimpl2(const T &eop, Clock &c) {
  c.start();
  int dummy = 0;
  double max;
  dso::MjdEpoch t(57752);
  while (t < dso::MjdEpoch(57754)) {
    dso::EopRecord ie;
    eop.interpolatev2(t, ie);
    max = (ie.xp() > ie.yp()) ? ie.xp() : ie.yp();
    max = (max > ie.dut()) ? max : ie.dut();
    max = (max > ie.lod()) ? max : ie.lod();
    max = (max > ie.dX()) ? max : ie.dX();
    t.add_seconds(30e0);
    /* flush cache */
    int i = (int)(max);
    dummy += flush_cache(i < 200 ? i : 200);
  }
  c.stop();
  c.add_sample();
  return dummy;
}

//int runimpl3(const dso::EopSeries &eop, Clock &c) {
//  c.start();
//  int dummy = 0;
//  double max;
//  dso::MjdEpoch t(57752);
//  while (t < dso::MjdEpoch(57754)) {
//    dso::EopRecord ie;
//    eop.interpolatev3(t, ie);
//    max = (ie.xp() > ie.yp()) ? ie.xp() : ie.yp();
//    max = (max > ie.dut()) ? max : ie.dut();
//    max = (max > ie.lod()) ? max : ie.lod();
//    max = (max > ie.dX()) ? max : ie.dX();
//    t.add_seconds(30e0);
//    /* flush cache */
//    int i = (int)(max);
//    dummy += flush_cache(i < 200 ? i : 200);
//  }
//  c.stop();
//  c.add_sample();
//  return dummy;
//}

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage %s [EOP C04 FILE]\n", argv[0]);
    return 1;
  }

  dso::MjdEpoch t1(57750);
  dso::MjdEpoch t2(57759);
  dso::EopSeries eop;

  if (dso::parse_iers_C04(argv[1], t1, t2, eop)) {
    fprintf(stderr, "ERROR Failed parsing eop file\n");
    return 1;
  }
  assert(eop.num_entries() == 9);
  dso::EopSeries2 eop2(eop);
  assert(eop2.num_entries() == 9);

  Clock c11("v11"), c21("v21"), c12("v12"),c22("v22");
  int dummy = 0;
  dummy += runimpl1(eop, c11);
  sleep(1);
  dummy += runimpl2(eop, c21);
  sleep(1);
  dummy += runimpl1(eop2, c12);
  sleep(1);
  dummy += runimpl2(eop2, c22);
  sleep(1);
  
  dummy += runimpl1(eop2, c12);
  sleep(1);
  dummy += runimpl1(eop, c11);
  sleep(1);
  dummy += runimpl2(eop, c21);
  sleep(1);
  dummy += runimpl2(eop2, c22);
  sleep(1);
  
  dummy += runimpl1(eop2, c12);
  sleep(1);
  dummy += runimpl2(eop, c21);
  sleep(1);
  dummy += runimpl2(eop2, c22);
  sleep(1);
  dummy += runimpl1(eop, c11);
  sleep(1);
  
  dummy += runimpl2(eop, c21);
  sleep(1);
  dummy += runimpl1(eop2, c12);
  sleep(1);
  dummy += runimpl2(eop2, c22);
  sleep(1);
  dummy += runimpl1(eop, c11);
  sleep(1);
  
  dummy += runimpl1(eop2, c12);
  sleep(1);
  dummy += runimpl2(eop2, c22);
  sleep(1);
  dummy += runimpl1(eop, c11);
  sleep(1);
  dummy += runimpl2(eop, c21);
  sleep(1);
  
  dummy += runimpl2(eop2, c22);
  sleep(1);
  dummy += runimpl1(eop, c11);
  sleep(1);
  dummy += runimpl1(eop2, c12);
  sleep(1);
  dummy += runimpl2(eop, c21);
  sleep(1);
  
  // choose faster clock
  Clock ref(c11);
  if (! ref.is_faster(c21) ) ref = c21;
  if (! ref.is_faster(c12) ) ref = c12;
  if (! ref.is_faster(c22) ) ref = c22;
  //if (! ref.is_faster(c3) ) ref = c3;
  c11.compare(ref);
  c21.compare(ref);
  c12.compare(ref);
  c22.compare(ref);
  //c3.compare(ref);
  printf("Average time for %s is %.2f\n", c11.mname, c11.maverage);
  printf("Average time for %s is %.2f\n", c12.mname, c12.maverage);
  printf("Average time for %s is %.2f\n", c21.mname, c21.maverage);
  printf("Average time for %s is %.2f\n", c22.mname, c22.maverage);

  printf("dummy = %d\n", dummy);
  return 0;
}
