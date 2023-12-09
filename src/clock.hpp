#include <chrono>
#include <cstring>
#include <stdint.h>
#include <cstdio>
using namespace std::chrono;

struct Clock {
  // uint_least64_t mruntime;
  double maverage;
  long msamples;
  time_point<high_resolution_clock> mstart;
  time_point<high_resolution_clock> mstop;
  char mname[24] = {'\0'};

  Clock(const char *n) { std::strcpy(mname, n); }

  Clock(const Clock &other) {
    // mruntime = other.mruntime;
    maverage = other.maverage;
    msamples = other.msamples;
    mstart = other.mstart;
    mstop = other.mstop;
    std::strcpy(mname, other.mname);
  }

  Clock(Clock &&other) {
    //mruntime = other.mruntime;
    maverage = other.maverage;
    msamples = other.msamples;
    mstart = other.mstart;
    mstop = other.mstop;
    std::strcpy(mname, other.mname);
  }

  Clock &operator=(const Clock &other) {
    if (this != &other) {
      //mruntime = other.mruntime;
      maverage = other.maverage;
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
    ++msamples;
    maverage = maverage * (msamples-1)/msamples + ((double)duration) / msamples;
    printf("%s average runtime=%.3f\n", mname, maverage);
  }
  
  bool is_faster(const Clock &other) const noexcept {
    return maverage < other.maverage;
  }
  
  void compare(const Clock &other) const noexcept {
    double ref = (double)(other.maverage);
    double y = (double)maverage;
    const double dx = 1e2 * (y-ref) / ref;
    printf("%s wrt %s : %+.1f %%\n", mname, other.mname, dx);
  }
};

