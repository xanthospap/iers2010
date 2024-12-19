target_sources(iers
  PRIVATE
    ${CMAKE_SOURCE_DIR}/src/lib/pole_tide/get_desai_ocp_deformation_coeffs.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/pole_tide/ocean_pole_tide_deformation.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/pole_tide/ocean_pole_tide_geopotential.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/pole_tide/parse_desai02_coeffs.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/pole_tide/pole_tide_deformation.cpp
)
