#include "kahan.hpp"
#include <cstdio>
#include <cassert>

using namespace dso;

int main() {
  KahanSum s;
  for (int i=0; i<1000; i++) s += 1e-3;

  assert( (double)s ==  1e0);

  return 0;
}
