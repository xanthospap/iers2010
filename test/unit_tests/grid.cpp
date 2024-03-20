#include "grid.hpp"
#include <cassert>
#include <cmath>

using namespace dso;

int main() {

  {
    GridAxis ga(2e0, 14.5e0, 2.5e0);
    ga.print();
    assert(ga.num_cells() == 5);
    double x = 2e0;
    while (x<=14.5e0) {
      auto x12 = ga.cell_nodes(ga.lower_bound(x));
      assert(x12.prev <= x && x12.next > x);
      int idx = ga.lower_bound(x);
      assert(ga.node(idx) <= x && x <= ga.node(idx+1));
      x += 1e-3;
    }
  }
  
  {
    GridAxis ga(-2e0, 19e0, 3e0);
    assert(ga.num_cells() == 7);
    double x = -2e0;
    while (x<=19e0) {
      auto x12 = ga.cell_nodes(ga.lower_bound(x));
      assert(x12.prev <= x && x12.next > x);
      int idx = ga.lower_bound(x);
      assert(ga.node(idx) <= x && x <= ga.node(idx+1));
      x += 1e-3;
    }
  }
  
  /* descending order */
  {
    GridAxis ga(19e0, -2e0, -3e0);
    assert(ga.num_cells() == 7);
    double x = 19e0;
    while (x>=-3e0) {
      auto x12 = ga.cell_nodes(ga.lower_bound(x));
      // assert(x12.prev <= x && x12.next > x);
      // printf("%.4f < %.4f < %.4f\n", x12.prev, x, x12.next);
      assert(x12.prev >= x && x12.next < x);
      int idx = ga.lower_bound(x);
      assert(ga.node(idx) >= x && x > ga.node(idx+1));
      x -= 1e-3;
    }
  }
  
  {
    GridAxis ga(-180e0, 180e0, 0.25e0);
    double x = -180e0;
    while (x<=180e0) {
      auto x12 = ga.cell_nodes(ga.lower_bound(x));
      assert(x12.prev <= x && x12.next > x);
      int idx = ga.lower_bound(x);
      assert(ga.node(idx) <= x && x <= ga.node(idx+1));
      x += 1e-3;
    }
  }
  
  {
    GridAxis ga(-M_PI, M_PI, 5e-5);
    double x = -M_PI;
    while (x<=M_PI) {
      auto x12 = ga.cell_nodes(ga.lower_bound(x));
      assert(x12.prev <= x && x12.next > x);
      int idx = ga.lower_bound(x);
      assert(ga.node(idx) <= x && x <= ga.node(idx+1));
      x += 1e-5;
    }
  }
  
  /* descending order */
  {
    GridAxis ga(M_PI, -M_PI, -5e-5);
    double x = M_PI;
    while (x>=-M_PI) {
      auto x12 = ga.cell_nodes(ga.lower_bound(x));
      // printf("%.9f >= %.9f > %.9f\n", x12.prev, x, x12.next);
      assert(x12.prev >= x && x12.next < x);
      int idx = ga.lower_bound(x);
      assert(ga.node(idx) >= x && x >= ga.node(idx+1));
      x -= 1e-5;
    }
  }

  return 0;
}
