#include "tropo.hpp"

void dso::gpt3::gpt3_grid::allocate(unsigned num_rows) {
  constexpr int num_cols = 5;
  p_grid = new double *[num_rows];
  for (unsigned i = 0; i < num_rows; i++)
    p_grid[i] = new double[num_cols];
  T_grid = new double *[num_rows];
  for (unsigned i = 0; i < num_rows; i++)
    T_grid[i] = new double[num_cols];
  Q_grid = new double *[num_rows];
  for (unsigned i = 0; i < num_rows; i++)
    Q_grid[i] = new double[num_cols];
  dT_grid = new double *[num_rows];
  for (unsigned i = 0; i < num_rows; i++)
    dT_grid[i] = new double[num_cols];
  u_grid = new double[num_rows];
  Hs_grid = new double[num_rows];
  ah_grid = new double *[num_rows];
  for (unsigned i = 0; i < num_rows; i++)
    ah_grid[i] = new double[num_cols];
  aw_grid = new double *[num_rows];
  for (unsigned i = 0; i < num_rows; i++)
    aw_grid[i] = new double[num_cols];
  la_grid = new double *[num_rows];
  for (unsigned i = 0; i < num_rows; i++)
    la_grid[i] = new double[num_cols];
  Tm_grid = new double *[num_rows];
  for (unsigned i = 0; i < num_rows; i++)
    Tm_grid[i] = new double[num_cols];
  Gn_h_grid = new double *[num_rows];
  for (unsigned i = 0; i < num_rows; i++)
    Gn_h_grid[i] = new double[num_cols];
  Ge_h_grid = new double *[num_rows];
  for (unsigned i = 0; i < num_rows; i++)
    Ge_h_grid[i] = new double[num_cols];
  Gn_w_grid = new double *[num_rows];
  for (unsigned i = 0; i < num_rows; i++)
    Gn_w_grid[i] = new double[num_cols];
  Ge_w_grid = new double *[num_rows];
  for (unsigned i = 0; i < num_rows; i++)
    Ge_w_grid[i] = new double[num_cols];
  size = num_rows;
}

void dso::gpt3::gpt3_grid::dealloc() noexcept {
  for (unsigned i = 0; i < size; i++)
    delete[] p_grid[i];
  delete[] p_grid;
  for (unsigned i = 0; i < size; i++)
    delete[] T_grid[i];
  delete[] T_grid;
  for (unsigned i = 0; i < size; i++)
    delete[] Q_grid[i];
  delete[] Q_grid;
  for (unsigned i = 0; i < size; i++)
    delete[] dT_grid[i];
  delete[] dT_grid;
  delete[] u_grid;
  delete[] Hs_grid;
  for (unsigned i = 0; i < size; i++)
    delete[] ah_grid[i];
  delete[] ah_grid;
  for (unsigned i = 0; i < size; i++)
    delete[] aw_grid[i];
  delete[] aw_grid;
  for (unsigned i = 0; i < size; i++)
    delete[] la_grid[i];
  delete[] la_grid;
  for (unsigned i = 0; i < size; i++)
    delete[] Tm_grid[i];
  delete[] Tm_grid;
  for (unsigned i = 0; i < size; i++)
    delete[] Gn_h_grid[i];
  delete[] Gn_h_grid;
  for (unsigned i = 0; i < size; i++)
    delete[] Ge_h_grid[i];
  delete[] Ge_h_grid;
  for (unsigned i = 0; i < size; i++)
    delete[] Gn_w_grid[i];
  delete[] Gn_w_grid;
  for (unsigned i = 0; i < size; i++)
    delete[] Ge_w_grid[i];
  delete[] Ge_w_grid;
  size = 0;
}
