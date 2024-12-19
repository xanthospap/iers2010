target_sources(iers
  PRIVATE
    ${CMAKE_SOURCE_DIR}/src/lib/gravity/forward_column_alf.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/gravity/modified_forward_column_alf.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/gravity/point_mass_acceleration.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/gravity/sh_cunningham_normalized_gradient.cpp
)
