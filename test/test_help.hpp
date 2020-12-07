#include <cassert>
#include <iostream>
#include <random>

inline 
double generate_random_double(double lower_bound, double upper_bound) noexcept {
  std::random_device rd;
  std::mt19937 e2(rd());
  std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
  return unif(e2);
}

// see
// https://stackoverflow.com/questions/17333/what-is-the-most-effective-way-for-float-and-double-comparison
template <typename TReal>
bool approxEqual(
    TReal a, TReal b,
    TReal tolerance = std::numeric_limits<TReal>::epsilon()) noexcept {
  TReal diff = std::abs(a - b);
  if (diff < tolerance)
    return true;
  if (diff < std::max(std::abs(a), std::abs(b)) * tolerance)
    return true;
  return false;
}

template <typename TReal>
bool approxEqual2(TReal a, TReal b, unsigned int interval_size = 1) noexcept {
  TReal min_a =
      a - (a - std::nextafter(a, std::numeric_limits<TReal>::lowest())) *
              interval_size;
  TReal max_a = a + (std::nextafter(a, std::numeric_limits<TReal>::max()) - a) *
                        interval_size;

  return min_a <= b && max_a >= b;
}
