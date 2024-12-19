target_sources(iers
  PRIVATE
    ${CMAKE_SOURCE_DIR}/src/lib/free_core_nutation/lambert.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/free_core_nutation/load_iers10_table52c.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/free_core_nutation/parse_lambert_coeffs.cpp
)
