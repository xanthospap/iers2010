#include "icgemio.hpp"
#include <stdexcept>

dso::Icgem::Icgem(const char *fn) : _filename(fn) {
  if (parse_header()) {
    throw std::runtime_error(
        "[ERROR] Failed to initialize input ICGEM file from " + _filename +
        "\n");
  }
}
