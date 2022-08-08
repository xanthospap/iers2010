#include "tropo.hpp"

using namespace dso::gpt3;

  int dso::Gpt3Grid::allocate(double resolution) {
    int new_rows = -1;
    if (resolution==1e0) {
      new_rows = gpt3_grid_attributes<Gpt3GridResolution::grid1x1>::num_lines;
    } else if (resolution==5e0) {
      new_rows = gpt3_grid_attributes<Gpt3GridResolution::grid5x5>::num_lines;
    } else {
      return 1;
    }

    if (new_rows != num_rows) {
      this->deallocate();
      memPool = new double[12*(5*new_rows) + 2*new_rows];
      offset = 5*new_rows;
      num_rows = new_rows;
    }

    return 0;
  }
