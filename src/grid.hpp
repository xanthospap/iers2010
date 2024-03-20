#ifndef __DSO_TOOLS_GRID_HPP__
#define __DSO_TOOLS_GRID_HPP__

#include <cstdio>
#include <cmath>
#include <stdexcept>

namespace dso {

/** @class GridAxis
 * A one-dimensional grid/axis that has 'nodes' and 'cells', as follows:
 *
 * node0   node1   node2                   nodeN
 * |       |       |                       |
 * V       V       V                       V
 * +-------+-------+-------+--...--+-------+
 * start           < step  >               stop
 * |cell0  | cell1 | cell2 |  ...  | celln |
 *
 */
class GridAxis {
  double _start;
  double _step;
  int _nodes;

public:
  struct CellNodes {
    double prev, next;
    CellNodes(double x1, double x2) noexcept : prev(x1), next(x2){};
  };

  /* Constructor 
   * Note that the value of to (passed in) may not be the value of the final 
   * node after construction. This is due to floating point precesion. The 
   * actual value of the final node will be:
   * from + _nodes * step
   * where
   * _nodes = (to - from) / step
   * */
  GridAxis(double from, double to, double step)
      : _start(from), _step(step), _nodes((to - from) / step) {
    if ((to - from < 0 && step > 0) || (to - from > 0 && step < 0)) {
      throw std::runtime_error(
          "[ERROR] Failed to construct GridAxis from given values!\n");
    }
  }

  /* number of nodes on axis, i.e. [from, to], both ends inclusive */
  int num_nodes() const noexcept { return _nodes + 1; }

  /* number of cells on axis */
  int num_cells() const noexcept { return num_nodes() - 1; }

  /* value of node with index idx */
  double node(int idx) const noexcept { return _start + idx * _step; }

  /* value of last node */
  double end() const noexcept { return _start + _nodes * _step; }

  /* If the axis is in ascending order (i.e. step>0), returns the index of the
   * node with the last value <= x, i.e. value of index+1 node will be > x.
   * int idx = lower_bound(123.456);
   * node(idx) <= 123.456 < node(idx+1)
   *
   * If the axis is in descending order (i.e. step<0), returns the index of
   * the node with the last value >= x, i.e. value of index+1 will be < x.
   * int idx = lower_bound(123.456);
   * node(idx) >= 123.456 > node(idx+1)
   *
   * Note that the limits of the instance are not considered here, and the 
   * value passed in could be outside the actual range of the instance (i.e. 
   * start/end nodes). It is the user's responsibility to check beforehand if 
   * the passed-in x value is indeed within the range of the instance.
   */
  int lower_bound(double x) const noexcept {
    return std::floor((x - _start) / _step);
  }

  int out_of_bounds(double x) const noexcept {
    return (x < _start || x > end());
  }

  /* Value of nodes of the cell. Nota that if the instance is in ascending 
   * order, then CellNodes.prev < CellNodes.next.
   * If the instance is in descending order, then 
   * CellNodes.prev > CellNodes.next.
   * 
   * Note that the limits of the instance are not considered here, and the 
   * value passed in could be outside the actual range of the instance (i.e. 
   * start/end nodes). It is the user's responsibility to check beforehand if 
   * the passed-in x value is indeed within the range of the instance.
   */
  CellNodes cell_nodes(int idx) noexcept {
    return CellNodes(node(idx), node(idx+1));
  }

  void print() const noexcept {
    printf("Cell from/to/step %.3f/%.3f/%.3f\n", _start * 1e0,
           _start * 1e0 + _nodes * _step, _step * 1e0);
  }
}; /* class GridAxis */

class TwoDimGrid {
  /* x is the first, leading dimension */
  GridAxis _x;
  /* y is the 'inner' dimension */
  GridAxis _y;
public:
  struct IndexPair {
    int xidx;
    int yidx;
    IndexPair(int i1, int i2) noexcept : xidx(i1), yidx(i2) {};
  }; /* IndexPair */

  TwoDimGrid(const GridAxis &x, const GridAxis &y) noexcept : _x(x), _y(y) {}
  
  /* number of nodes on grid */
  int num_nodes() const noexcept { return _x.num_nodes() * _y.num_nodes(); }
  int num_x_nodes() const noexcept { return _x.num_nodes();}
  int num_y_nodes() const noexcept { return _y.num_nodes();}

  GridAxis &xaxis() noexcept { return _x; }
  const GridAxis &xaxis() const noexcept { return _x; }
  GridAxis &yaxis() noexcept { return _y; }
  const GridAxis &yaxis() const noexcept { return _y; }

}; /* class TwoDimGrid */
} /* namespace dso */

#endif
