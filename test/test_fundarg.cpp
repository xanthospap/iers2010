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

  // results as recorded in the FORTRAN source file for the test case
  double fargs_ref[] = {2.291187512612069099e0, 6.212931111003726414e0,
                        3.658025792050572989e0, 4.554139562402433228e0,
                        -0.5167379217231804489e0};

  // discrepancies between the FORTRAN implementation and fargs_ref; should be
  // the same as C++
  double fortran_diffs[] = {0.77e-13, 0.97e-11, 0.78e-13, 0.96e-11, 0.23e-15};

  iers2010::fundarg(0.07995893223819302e0, fargs);
  for (int i = 0; i < 5; ++i) {
    printf("\nargs[%1d] = %12.6e radians", i,
           std::abs(fargs[i] - fargs_ref[i]));
    assert(std::abs(fargs[i] - fargs_ref[i]) < _alg_accuracy_);
    assert(std::abs(fargs[i] - fargs_ref[i]) < std::abs(fortran_diffs[i]));
  }
  return 0;
}
