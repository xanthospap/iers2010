#include <chrono>
#include <cstring>
#include <stdint.h>
#include <cstdio>

struct Clock {
  using clock_t = std::chrono::high_resolution_clock;
  using duration_t = clock_t::duration;
  
  clock_t::time_point start_time;
  clock_t::duration total_duration{0};
  int msamples{0};
  char mname[24] = {'\0'};

  duration_t get_elapsed_time() const { return (clock_t::now() - start_time); }
  void start() { start_time = clock_t::now(); }
  void stop() {
    total_duration += get_elapsed_time();
    ++msamples;
  }
  template <typename T = std::chrono::microseconds> auto get_total_duration() const {
    return std::chrono::duration_cast<T>(total_duration);
  }
  template <typename T = std::chrono::microseconds>
  double average_running_time() const {
    return get_total_duration<T>().count() / msamples;
  }

  Clock(const char *n) { std::strcpy(mname, n); }

  Clock(const Clock &other) {
    start_time = other.start_time;
    total_duration = other.total_duration;
    msamples = other.msamples;
    std::strcpy(mname, other.mname);
  }

  Clock(Clock &&other) {
    start_time = other.start_time;
    total_duration = other.total_duration;
    msamples = other.msamples;
    std::strcpy(mname, other.mname);
  }

  Clock &operator=(const Clock &other) {
    if (this != &other) {
      start_time = other.start_time;
      total_duration = other.total_duration;
      msamples = other.msamples;
      std::strcpy(mname, other.mname);
    }
    return *this;
  }

  const char *name() const {return mname;}

  bool is_faster(const Clock &other) const noexcept {
    return average_running_time() < other.average_running_time();
  }
  
  void compare(const Clock &other) const noexcept {
    double ref = other.average_running_time();
    double y = average_running_time();
    const double dx = 1e2 * (y-ref) / ref;
    printf("%s wrt %s : %+.1f %%\n", mname, other.mname, dx);
  }
};

