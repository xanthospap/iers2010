#include <fstream>
#include <cstdio>
#include <cstring>
#include "gpt3.hpp"

int dso::gpt3::parse_gpt3_5_grid(const char *gridfn, dso::gpt3::gpt3_5_grid *grid) noexcept {
  std::ifstream grd(gridfn);
  if (!grd) {
    fprintf(stderr,
            "[ERROR] Failed to open gpt3_5 grid file %s (traceback: %s)\n",
            gridfn, __func__);
    return 1;
  }

  // read first line, that is the header
  char header[512];
  grd.getline(header, 512);
  if (!grd.good()) {
    fprintf(stderr,
            "[ERROR] Failed extracting first/header line from gpt3_5 grid file "
            "%s (traceback: %s)\n",
            gridfn, __func__);
    return 2;
  }

  // expected header
  const char *chdr =
      "%  lat    lon   p:a0    A1   B1   A2   B2  T:a0   A1    B1   A2   B2  "
      "Q:a0   A1    B1    A2    B2    dT:a0   A1   B1   A2   B2  undu      Hs  "
      "  a_h:a0     A1       B1        A2      B2     a_w:a0     A1       B1   "
      "    A2       B2    lambda:a0  A1     B1      A2     B2    Tm:a0   A1   "
      "B1   A2   B2    Gn_h:a0   A1    B1     A2     B2   Ge_h:a0   A1    B1   "
      "  A2     B2   Gn_w:a0   A1    B1     A2     B2   Ge_w:a0   A1    B1     "
      "A2     B2";
  int chdr_len = std::strlen(chdr);

    // validate header
  if (std::strncmp(chdr, header, chdr_len)) {
    fprintf(stderr,
            "[ERROR] Failed to valildate header for gpt3_5 grid file %s "
            "(traceback: %s)\n",
            gridfn, __func__);
    return 3;
  }

    // read every data row ...
    double dm1, dm2;
  int row = 0;
  while (!grd.eof()) {
    grd >> dm1 >> dm2 >> grid->p_grid[row][0] >> grid->p_grid[row][1] >> grid->p_grid[row][2] >>
        grid->p_grid[row][3] >> grid->p_grid[row][4] >> grid->T_grid[row][0] >>
        grid->T_grid[row][1] >> grid->T_grid[row][2] >> grid->T_grid[row][3] >>
        grid->T_grid[row][4] >> grid->Q_grid[row][0] >> grid->Q_grid[row][1] >>
        grid->Q_grid[row][2] >> grid->Q_grid[row][3] >> grid->Q_grid[row][4] >>
        grid->dT_grid[row][0] >> grid->dT_grid[row][1] >> grid->dT_grid[row][2] >>
        grid->dT_grid[row][3] >> grid->dT_grid[row][4] >> grid->u_grid[row] >>
        grid->Hs_grid[row] >> grid->ah_grid[row][0] >> grid->ah_grid[row][1] >>
        grid->ah_grid[row][2] >> grid->ah_grid[row][3] >> grid->ah_grid[row][4] >>
        grid->aw_grid[row][0] >> grid->aw_grid[row][1] >> grid->aw_grid[row][2] >>
        grid->aw_grid[row][3] >> grid->aw_grid[row][4] >> grid->la_grid[row][0] >>
        grid->la_grid[row][1] >> grid->la_grid[row][2] >> grid->la_grid[row][3] >>
        grid->la_grid[row][4] >> grid->Tm_grid[row][0] >> grid->Tm_grid[row][1] >>
        grid->Tm_grid[row][2] >> grid->Tm_grid[row][3] >> grid->Tm_grid[row][4] >>
        grid->Gn_h_grid[row][0] >> grid->Gn_h_grid[row][1] >>
        grid->Gn_h_grid[row][2] >> grid->Gn_h_grid[row][3] >>
        grid->Gn_h_grid[row][4] >> grid->Ge_h_grid[row][0] >>
        grid->Ge_h_grid[row][1] >> grid->Ge_h_grid[row][2] >>
        grid->Ge_h_grid[row][3] >> grid->Ge_h_grid[row][4] >>
        grid->Gn_w_grid[row][0] >> grid->Gn_w_grid[row][1] >>
        grid->Gn_w_grid[row][2] >> grid->Gn_w_grid[row][3] >>
        grid->Gn_w_grid[row][4] >> grid->Ge_w_grid[row][0] >>
        grid->Ge_w_grid[row][1] >> grid->Ge_w_grid[row][2] >>
        grid->Ge_w_grid[row][3] >> grid->Ge_w_grid[row][4];
    ++row;
  }

  if (row != GPT3_5_GRID_LINES) {
    fprintf(
        stderr,
        "[ERROR] Parsed %d/%d lines from gpt3_5 grid file %s (traceback: %s)\n",
        row, GPT3_5_GRID_LINES, gridfn, __func__);
    return 10;
  }

  // handle units ...
  for (int i=0; i<GPT3_5_GRID_LINES; i++) for (int j=0; j<5; j++) grid->Q_grid[i][j] /= 1e3;
  for (int i=0; i<GPT3_5_GRID_LINES; i++) for (int j=0; j<5; j++) grid->dT_grid[i][j] /= 1e3;
  for (int i=0; i<GPT3_5_GRID_LINES; i++) for (int j=0; j<5; j++) grid->ah_grid[i][j] /= 1e3;
  for (int i=0; i<GPT3_5_GRID_LINES; i++) for (int j=0; j<5; j++) grid->aw_grid[i][j] /= 1e3;
  for (int i=0; i<GPT3_5_GRID_LINES; i++) for (int j=0; j<5; j++) grid->Gn_h_grid[i][j] /= 1e5;
  for (int i=0; i<GPT3_5_GRID_LINES; i++) for (int j=0; j<5; j++) grid->Ge_h_grid[i][j] /= 1e5;
  for (int i=0; i<GPT3_5_GRID_LINES; i++) for (int j=0; j<5; j++) grid->Gn_w_grid[i][j] /= 1e5;
  for (int i=0; i<GPT3_5_GRID_LINES; i++) for (int j=0; j<5; j++) grid->Ge_w_grid[i][j] /= 1e5;

  return 0;
}