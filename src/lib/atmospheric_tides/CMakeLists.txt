target_sources(iers
  PRIVATE
    ${CMAKE_SOURCE_DIR}/src/lib/atmospheric_tides/atm_tides_append_wave.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/atmospheric_tides/atm_tides_geopotential.cpp
)
