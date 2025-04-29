target_sources(iers
  PRIVATE
    ${CMAKE_SOURCE_DIR}/src/lib/earth_rotation/gcrs2itrs_implementation.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/earth_rotation/gcrs2itrs_quaternion.cpp
)
