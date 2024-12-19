#ifndef __DSO_TOOLS_GRID_HPP__
#define __DSO_TOOLS_GRID_HPP__

#include <cmath>
#include <cstdio>

namespace dso {
enum class GridAxis { X, Y };

template <GridAxis LeadingAxis> class TwoDimGrid {
public:
  struct Node {
    double x, y;
  };

private:
  int mxpts, mypts;
  double mxstart, mxstep;
  double mystart, mystep;

  int xindex(double xpts) const noexcept {
    return std::floor((xpts - mxstart) / mxstep);
  }

  int yindex(double ypts) const noexcept {
    return std::floor((ypts - mystart) / mystep);
  }

  long xyindex(double xpts, double ypts) const noexcept {
    auto xidx = xindex(xpts);
    auto yidx = yindex(ypts);
    if constexpr (LeadingAxis == GridAxis::X) {
      return xidx * mypts + yidx;
    } else {
      return yidx * mxpts + xidx;
    }
  }

  bool out_of_range_x(double x) const noexcept {
    return !(x >= mxstart && x < xstop());
  }
  bool out_of_range_y(double y) const noexcept {
    return !(y >= mystart && y < ystop());
  }

  void surrounding_nodes(long idx_bl, long &idx_br, long &idx_tl,
                         long &idx_tr) const noexcept {
    if constexpr (LeadingAxis == GridAxis::X) {
      idx_br = idx_bl + mypts;
      idx_tl = idx_bl + 1;
      idx_tr = idx_tl + 1;
    } else {
      idx_br = idx_bl + 1;
      idx_tl = idx_bl + mxpts;
      idx_tr = idx_tl + 1;
    }
  }

public:
  TwoDimGrid() noexcept
      : mxpts(0), mypts(0), mxstart(0), mxstep(0), mystart(0), mystep(0) {};
  TwoDimGrid(double xstart, double xstep, int xpts, double ystart, double ystep,
             int ypts)
      : mxpts(xpts), mypts(ypts), mxstart(xstart), mxstep(xstep),
        mystart(ystart), mystep(ystep) {};

  double xstep() const noexcept { return mxstep; }
  double ystep() const noexcept { return mystep; }
  double xstop() const noexcept { return mxstart + mxpts * mxstep; }
  double ystop() const noexcept { return mystart + mypts * mystep; }

  long line(double x, double y) const noexcept {
    if (out_of_range_x(x) || out_of_range_y(y)) {
      fprintf(stderr,
              "[ERROR] Given point is out of range of grid; (%.3f, %.3f) for "
              "[%.2f/%.2f/%.2f/%.2f] "
              "(traceback: %s)\n",
              x, y, mxstart, xstop(), mystart, ystop(), __func__);
      return -1;
    }
    return xyindex(x, y);
  }

  Node xyindex2node(long index) const noexcept {
    if constexpr (LeadingAxis == GridAxis::X) {
      long xi = index / mypts;
      long yi = index - xi * mxpts;
      return Node{xi * mxstep + mxstart, yi * mystep + mystart};
    } else {
      int yi = index / mxpts;
      int xi = index - yi * mxpts;
      return Node{xi * mxstep + mxstart, yi * mystep + mystart};
    }
  }

  void surrounding_nodes(double x, double y, long &idx_bl, long &idx_br,
                         long &idx_tl, long &idx_tr) const noexcept {
    idx_bl = line(x, y);
    if (idx_bl < 0)
      idx_bl = idx_br = idx_tl = idx_tr = -1;
    surrounding_nodes(idx_bl, idx_br, idx_tl, idx_tr);
  }

}; /* class TwoDimGrid */
} /* namespace dso */

#endif
