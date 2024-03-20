#include "datetime/calendar.hpp"
#include "eigen3/Eigen/Eigen"
#include "gravity.hpp"
#include "icgemio.hpp"
#include <charconv>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <vector>

struct BmOrbit {
  dso::MjdEpoch epoch;
  Eigen::Matrix<double, 3, 1> xyz;
  Eigen::Matrix<double, 3, 1> vxyz;
  Eigen::Matrix<double, 3, 1> axyz;
  BmOrbit(const dso::MjdEpoch &t, double *rva) : epoch(t) {
    xyz << rva[0], rva[1], rva[2];
    vxyz << rva[3], rva[4], rva[5];
    axyz << rva[6], rva[7], rva[8];
  }
};

struct BmAcceleration {
  dso::MjdEpoch epoch;
  Eigen::Matrix<double, 3, 1> axyz;
  BmAcceleration(const dso::MjdEpoch &t, double ax, double ay, double az)
      : epoch(t) {
    axyz << ax, ay, az;
  }
};

const char *skipws(const char *line) noexcept {
  while (line && *line == ' ')
    ++line;
  return line;
}

std::vector<BmAcceleration> parse_acceleration(const char *fn) {
  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR Failed opening file %s\n", fn);
    throw std::runtime_error("Failed opening input file " + std::string(fn) +
                             "\n");
  }

  char line[256];
  for (int i = 0; i < 5; i++)
    fin.getline(line, 256);
  std::vector<BmAcceleration> vec;

  int error = 0;
  while (fin.getline(line, 256) && (!error)) {
    double td[4];
    int sz = std::strlen(line);
    const char *str = line;
    for (int i = 0; i < 4; i++) {
      auto res = std::from_chars(skipws(str), line + sz, td[i]);
      if (res.ec != std::errc{})
        ++error;
      str = res.ptr;
    }
    if (!error) {
      int imjd = (int)td[0];
      double fsec = (td[0] - imjd) * 86400e0;
      const dso::MjdEpoch t(imjd, dso::FractionalSeconds{fsec});
      vec.emplace_back(t, td[1], td[2], td[3]);
    }
  }

  if (error || (!fin.eof())) {
    fprintf(stderr, "ERROR Failed parsing input file %s\n", fn);
    throw std::runtime_error("Failed parsing input file " + std::string(fn) +
                             "\n");
  }

  return vec;
}

std::vector<BmOrbit> parse_orbit(const char *fn) {
  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR Failed opening file %s\n", fn);
    throw std::runtime_error("Failed opening input file " + std::string(fn) +
                             "\n");
  }

  char line[512];
  for (int i = 0; i < 5; i++)
    fin.getline(line, 512);
  std::vector<BmOrbit> vec;

  int error = 0;
  while (fin.getline(line, 512) && (!error)) {
    double td[10];
    int sz = std::strlen(line);
    const char *str = line;
    for (int i = 0; i < 10; i++) {
      auto res = std::from_chars(skipws(str), line + sz, td[i]);
      if (res.ec != std::errc{})
        ++error;
      str = res.ptr;
    }
    if (!error) {
      int imjd = (int)td[0];
      double fsec = (td[0] - imjd) * 86400e0;
      const dso::MjdEpoch t(imjd, dso::FractionalSeconds{fsec});
      vec.emplace_back(t, &td[1]);
    }
  }

  if (error || (!fin.eof())) {
    fprintf(stderr, "ERROR Failed parsing input file %s\n", fn);
    throw std::runtime_error("Failed parsing input file " + std::string(fn) +
                             "\n");
  }

  return vec;
}

constexpr const int DEGREE = 180;
constexpr const int ORDER = 180;
int main(int argc, char *argv[]) {
  if (argc != 4) {
    fprintf(
        stderr,
        "Usage: %s [00orbit_itrf.txt] [GEOPOTENTIAL] [02gravityfield_itrf]\n",
        argv[0]);
    return 1;
  }

  /* read orbit from input file */
  const auto orbvec = parse_orbit(argv[1]);
  printf("#Note: read %d data sets from input file\n", (int)orbvec.size());

  /* read accleration from input file */
  const auto accvec = parse_acceleration(argv[3]);
  printf("#Note: read %d data sets from input file\n", (int)accvec.size());

  /* read gravity model into a StokesCoeffs instance */
  dso::Icgem icgem(argv[2]);
  dso::StokesCoeffs sc(DEGREE, ORDER, 0e0, 0e0);
  dso::Datetime<dso::nanoseconds> t(
      dso::from_mjdepoch<dso::nanoseconds>(orbvec[0].epoch));
  if (icgem.parse_data(DEGREE, ORDER, t, sc)) {
    fprintf(stderr, "ERROR Failed reading gravity model!\n");
    return 1;
  }
  // printf("C(180,180)=%.15e\n", sc.C(180,180));
  // assert(-0.22147242139e-09 == sc.C(180, 180));
  assert(-0.69378736267e-09 == sc.C(150, 150));
  // return 10;

  /* allocate scratch space for computations */
  dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> W(DEGREE + 3,
                                                                    DEGREE + 3);
  dso::CoeffMatrix2D<dso::MatrixStorageType::LwTriangularColWise> M(DEGREE + 3,
                                                                    DEGREE + 3);

  /* compare results epoch by epoch */
  Eigen::Matrix<double, 3, 1> a;
  Eigen::Matrix<double, 3, 3> g;
  auto acc = accvec.begin();
  for (const auto &in : orbvec) {
    /* compute acceleration for given epoch/position */
    if (dso::sh2gradient_cunningham(sc, in.xyz, a, g, DEGREE, ORDER, -1, -1, &W,
                                    &M)) {
      fprintf(stderr, "ERROR Failed computing acceleration/gradient\n");
      return 1;
    }
    /* get COSTG result */
    if (acc->epoch != in.epoch) {
      fprintf(stderr, "ERROR Faile to match epochs in input files\n");
      return 1;
    }
    printf("%d %.9f %.15e %.15e %.15e %.15e %.15e %.15e\n", in.epoch.imjd(),
           in.epoch.seconds(), acc->axyz(0), acc->axyz(1), acc->axyz(2), a(0),
           a(1), a(2));
    ++acc;
  }

  return 0;
}
