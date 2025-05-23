#include "aod1b_data_stream.hpp"
#include "datetime/datetime_write.hpp"
#include <cstdio>
#include <datetime/fractional_types.hpp>
#include <iers/stokes_coefficients.hpp>
#include <stdexcept>
#include <random>

int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s [AOD1B_FILE] [AOD1B_DIR]\n", argv[0]);
    return 1;
  }

  std::uniform_real_distribution<double> unif(0., 3600e0);
  std::default_random_engine re;

  dso::StokesCoeffs s(80, 80);

  [[maybe_unused]] char buf[64];

  dso::Aod1bDataStream<dso::AOD1BCoefficientType::GLO> *dap = nullptr;

    try {
      dap = new dso::Aod1bDataStream<dso::AOD1BCoefficientType::GLO>(argv[1], argv[2]);
      dap->initialize();
    } catch (std::exception &e) {
      fprintf(stderr, e.what());
      std::string err_msg =
          "[ERROR] Failed to construct an Aod1bDataStream<GLO> stream from "
          "config parameters " +
          std::string(argv[1]) + ", " + std::string(argv[2]) + "\n";
      throw std::runtime_error(err_msg);
    }

    auto t = dap->stream().first_epoch().add_seconds(dso::nanoseconds(unif(re)*1e9));
    const auto maxt = dap->stream().first_epoch().add_seconds(dso::seconds(86400+86400/2));
    double advanced = 0;
    while (t < maxt) {
      if (dap->coefficients_at(t, s)) return 100;
      t >= dap->
      // printf("got dealiasing for %s\n",dso::to_char<dso::YMDFormat::YYYYMMDD, dso::HMSFormat::HHMMSSF, dso::nanoseconds>(t, buf));
      double nsec = unif(re);
      if (advanced>3600.) {
        while (nsec<30.) nsec = unif(re);
        advanced = 0.;
        nsec *= -1e0;
      } else {
        advanced += nsec;
      }
      t.add_seconds_inplace(dso::nanoseconds(nsec*1e9));
    }

  return 0;
}
