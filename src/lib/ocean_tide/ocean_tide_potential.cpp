#include "iau.hpp"
#include "ocean_tide.hpp"
#include "eigen3/Eigen/Dense"
#include <iostream>

int dso::OceanTide::stokes_coeffs_impl(const dso::MjdEpoch& mjdtt,
    const dso::MjdEpoch& mjdut1,
    const double* const delaunay_args) noexcept {
  /* nullify geopotential coeffs */
  mcs.clear();

  /* compute GMST using IAU 2006/2000A [rad] */
  const double gmst = dso::gmst(mjdtt, mjdut1);

  /* compute six-vector of multipliers ni from Delaunay vars */
  double __dargs[6];
  const double* __restrict__ f = dso::delaunay2doodson(delaunay_args, gmst, __dargs);

  /* iterate through individual constituents */
  for (const auto& wave : atlas().waves()) {
    /* compute angle: θ(f) = Σ(i=1,6) n(i)*β(i) */
    const double arg = wave.wave().doodson().argument(f);
    const double carg = std::cos(arg);
    const double sarg = std::sin(arg);
    mcs.Cnm() += wave.stokes_cos().Cnm() * carg + wave.stokes_sin().Cnm() * sarg;
    mcs.Snm() += wave.stokes_cos().Snm() * carg + wave.stokes_sin().Snm() * sarg;
  }

  return 0;
}

int dso::OceanTide::stokes_coeffs_with_admittance_impl(const dso::MjdEpoch& mjdtt,
    const dso::MjdEpoch& mjdut1,
    const double* const delaunay_args) noexcept {
  /* nullify geopotential coeffs */
  mcs.clear();

  /* compute GMST using IAU 2006/2000A [rad] */
  const double gmst = dso::gmst(mjdtt, mjdut1);

  /* compute six-vector of multipliers ni from Delaunay vars */
  //double __dargs[6];
  //dso::delaunay2doodson(delaunay_args, gmst, __dargs);
  //Eigen::Matrix<double, 6, 1> theta = Eigen::Map<Eigen::Matrix<double, 6, 1>>(__dargs);

  /* compute theta angles for all waves (fx1) */
  //const Eigen::VectorXd thetaf = doodson_matrix().cast<double>() * theta;
  
  double __dargs[6];
  const double* __restrict__ f = dso::delaunay2doodson(delaunay_args, gmst, __dargs);
  Eigen::VectorXd thetaf = Eigen::VectorXd::Zero(doodson_matrix().rows());
  int _dint[6];
  for (int i=0; i<doodson_matrix().rows(); i++) {
    for (int j=0; j<6; j++) _dint[j] = doodson_matrix()(i,j);
    dso::DoodsonConstituent ds(&_dint[0]);
    thetaf(i) = ds.argument(f);
  }

  /* factors for cos and sin (kx1) */
  //std::cout<<"admittance:\n"<<admittance_matrix()<<"\n";
  const Eigen::VectorXd factorCos = admittance_matrix() * thetaf.unaryExpr([](double x){return std::cos(x);});
  const Eigen::VectorXd factorSin = admittance_matrix() * thetaf.unaryExpr([](double x){return std::sin(x);});
  //std::cout << "factorCos=\n" << factorCos << "\n";

  /* iterate through individual constituents */
  const double* __restrict__ fcos = factorCos.data();
  const double* __restrict__ fsin = factorSin.data();
  for (const auto& wave : atlas().waves()) {
    mcs.Cnm() += wave.stokes_cos().Cnm() * (*fcos) + wave.stokes_sin().Cnm() * (*fsin);
    mcs.Snm() += wave.stokes_cos().Snm() * (*fcos) + wave.stokes_sin().Snm() * (*fsin);
    //printf("fcos = %.12e fsin = %.12e\n", *fcos, *fsin);
    ++fcos;
    ++fsin;
  }

  return 0;
}
