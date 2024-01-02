#include "aod1b_data_stream.hpp"
#include "datetime/datetime_write.hpp"
#include <cstdio>

constexpr const dso::nanoseconds _01_hour_nanosec (60L     * 60L * 1'000'000'000L);
constexpr const dso::nanoseconds _03_hour_nanosec (3 * 60L * 60L * 1'000'000'000L);
constexpr const dso::nanoseconds _07_hour_nanosec (7 * 60L * 60L * 1'000'000'000L);
constexpr const dso::nanoseconds _23_hour_nanosec (23 * 60L* 60L * 1'000'000'000L);

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s [AOD1B]\n", argv[0]);
    return 1;
  }

  dso::Aod1bIn aod(argv[1]);
  dso::Aod1bDataStream<dso::AOD1BCoefficientType::ATM> stream(aod);
  dso::StokesCoeffs cs(120,120,0e0,0e0);

  if (stream.initialize()) {
    return 1;
  }
  auto t = aod.first_epoch();
  while (t < aod.last_epoch()) {
    if (stream.coefficients_at(t,cs)) return 8;
    t.add_seconds(dso::seconds(1));
  }

  return 0;
}
