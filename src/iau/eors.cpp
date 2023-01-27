#include "iau.hpp"

double iers2010::sofa::eors(const Eigen::Matrix<double,3,3> &rnpb, double s) noexcept {
  // Evaluate Wallace & Capitaine (2006) expression (16).
  const double x = rnpb(2,0); // [2][0]
  const double ax = x / (1e0 + rnpb(2,2)/*[2][2]*/);
  const double xs = 1e0 - ax * x;
  const double ys = -ax * rnpb(2,1); // [2][1]
  const double zs = -x;
  const double p = rnpb(0,0) /*[0][0]*/ * xs + rnpb(0,1) /*[0][1]*/ * ys +
                   rnpb(0,2) /*[0][2]*/ * zs;
  const double q = rnpb(1,0) /*[1][0]*/ * xs + rnpb(1,1) /*[1][1]*/ * ys +
                   rnpb(1,2) /*[1][2]*/ * zs;
  const double eo = ((p != 0) || (q != 0)) ? s - std::atan2(q, p) : s;

  return eo;
}
