target_sources(iers
  PRIVATE
    ${CMAKE_SOURCE_DIR}/src/lib/earth_rotation/gcrs_to_itrs.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/earth_rotation/icrs2itrs_implementation.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/earth_rotation/itrs2gcrs_quaternion.cpp
)
