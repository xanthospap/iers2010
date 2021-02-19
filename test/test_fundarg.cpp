#include "iers2010.hpp"
#include "test_help.hpp"
#include <cassert>

constexpr double _alg_accuracy_ = 1e-11;

int main() {

  // testing fundarg
  std::cout << "----------------------------------------\n";
  std::cout << "> fundarg\n";
  std::cout << "----------------------------------------\n";

  double fargs[5];
  double fargs_ref[] = {2.291187512612069099e0, 6.212931111003726414e0,
                        3.658025792050572989e0, 4.554139562402433228e0,
                        -0.5167379217231804489e0};

  iers2010::fundarg(0.07995893223819302e0, fargs);
  for (int i = 0; i < 5; ++i) {
#ifdef STRICT_TEST
    assert(approxEqual(fargs[i], fargs_ref[i]));
#else
    printf("\nargs[%1d] = %12.6e radians", i,
           std::abs(fargs[i] - fargs_ref[i]));
    assert(std::abs(fargs[i] - fargs_ref[i]) < _alg_accuracy_);
#endif
  }
  return 0;
}
