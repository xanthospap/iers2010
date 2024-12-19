target_sources(iers
  PRIVATE
    ${CMAKE_SOURCE_DIR}/src/lib/stokes_coefficients.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/tidal_constituent.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/tide_atlas.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/tide_atlas_groops.cpp
)
