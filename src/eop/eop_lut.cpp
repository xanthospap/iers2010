#include "eop.hpp"
  
dso::EopLookUpTable::~EopLookUpTable() noexcept = default;

void dso::EopLookUpTable::clear() noexcept {
  t.clear();
  xp.clear();
  yp.clear();
  dut1.clear();
  dX.clear();
  dY.clear();
  lod.clear();
  xrt.clear();
  yrt.clear();
}

void dso::EopLookUpTable::reserve(int sz) noexcept {
  t.reserve(sz);
  xp.reserve(sz);
  yp.reserve(sz);
  dut1.reserve(sz);
  dX.reserve(sz);
  dY.reserve(sz);
  lod.reserve(sz);
  xrt.reserve(sz);
  yrt.reserve(sz);
}

void dso::EopLookUpTable::push_back(const dso::EopRecord& rec) noexcept {
  t.push_back(rec.mjd);
  xp.push_back(rec.xp);
  yp.push_back(rec.yp);
  dut1.push_back(rec.dut);
  dX.push_back(rec.dx);
  dY.push_back(rec.dy);
  lod.push_back(rec.lod);
  xrt.push_back(rec.xrt);
  yrt.push_back(rec.yrt);
}
