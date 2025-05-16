target_sources(iers
  PRIVATE
    ${CMAKE_SOURCE_DIR}/src/lib/earth_rotation/c2i06a.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/earth_rotation/c2i06abz.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/earth_rotation/c2i06apq.cpp
)
