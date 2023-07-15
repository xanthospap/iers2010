#include "eop.hpp"
  
dso::EopLookUpTable::~EopLookUpTable() noexcept = default;

void dso::EopLookUpTable::clear() noexcept {
  tvec.clear();
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
  tvec.reserve(sz);
  xp.reserve(sz);
  yp.reserve(sz);
  dut1.reserve(sz);
  dX.reserve(sz);
  dY.reserve(sz);
  lod.reserve(sz);
  xrt.reserve(sz);
  yrt.reserve(sz);
}

void dso::EopLookUpTable::push_back(const dso::EopRecord &rec) noexcept {
  /* get an iterator to the upper bound, i.e. it > t */
  auto it = this->upper_bound(rec.mjd);

  if (it == tvec.end()) {
    /* simplest case, new element to be added at end */
    tvec.push_back(rec.mjd);
    xp.push_back(rec.xp);
    yp.push_back(rec.yp);
    dut1.push_back(rec.dut);
    dX.push_back(rec.dx);
    dY.push_back(rec.dy);
    lod.push_back(rec.lod);
    xrt.push_back(rec.xrt);
    yrt.push_back(rec.yrt);
  } else {
    /* inserts value(s) just before before the iterator */
    tvec.insert(it, rec.mjd);
    /* get the index of the iterator */
    const int index = std::distance(tvec.begin(), it);
    xp.insert(xp.begin() + index, rec.xp);
    yp.insert(yp.begin() + index, rec.yp);
    dut1.insert(dut1.begin() + index, rec.dut);
    dX.insert(dX.begin() + index, rec.dx);
    dY.insert(dY.begin() + index, rec.dy);
    lod.insert(lod.begin() + index, rec.lod);
    xrt.insert(xrt.begin() + index, rec.xrt);
    yrt.insert(yrt.begin() + index, rec.yrt);
  }

  return;
}
