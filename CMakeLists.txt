cmake_minimum_required(VERSION 3.8)

project(
  iers
  VERSION 1.1.0
  DESCRIPTION "Implementation of IERS2010 (upd.18) standards"
  LANGUAGES CXX
)

# Ensure required libraries are available
find_package(Eigen3   REQUIRED)
find_package(geodesy  REQUIRED)
find_package(datetime REQUIRED)
find_library(CSPICE_LIB cspice REQUIRED)

# Pass the library dependencies to subdirectories
set(PROJECT_DEPENDENCIES Eigen3::Eigen geodesy datetime ${CSPICE_LIB})

# Define an option for building tests (defaults to ON)
option(BUILD_TESTING "Enable building of tests" OFF)

# Define the option to build costg
option(BUILD_COSTG "Build costg binaries (requires iers library)" OFF)

# compiler flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED On)
set(CMAKE_CXX_EXTENSIONS Off)
add_compile_options(-Wall 
  -Wextra 
  -Werror 
  -pedantic 
  -W 
  -Wshadow 
  $<$<CONFIG:Release>:-O2>
  $<$<CONFIG:Release>:-march=native>
  $<$<CONFIG:Debug>:-g>
  $<$<CONFIG:Debug>:-pg> 
  $<$<CONFIG:Debug>:-Wdisabled-optimization>
  )
add_compile_definitions(
  $<$<CONFIG:Debug>:DEBUG>
)

# the library and includes
add_library(iers)
target_include_directories(iers
PUBLIC
  # Public headers during build time
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
  # Public headers during install time
  $<INSTALL_INTERFACE:include/iers>
PRIVATE
  # Private/internal headers only during build time (src/core)
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  # Internal/private headers during install time (from include/iers/core)
  $<INSTALL_INTERFACE:include/iers/core>
)

# library source code
add_subdirectory(src/lib)
add_subdirectory(src/lib/aod1b)
add_subdirectory(src/lib/atmospheric_tides)
add_subdirectory(src/lib/doodson)
add_subdirectory(src/lib/earth_rotation)
add_subdirectory(src/lib/eop)
add_subdirectory(src/lib/free_core_nutation)
add_subdirectory(src/lib/fundarg)
add_subdirectory(src/lib/gravity)
add_subdirectory(src/lib/iau)
add_subdirectory(src/lib/icgem)
add_subdirectory(src/lib/ocean_tide)
add_subdirectory(src/lib/planets)
add_subdirectory(src/lib/pole_tide)
add_subdirectory(src/lib/relativity)
add_subdirectory(src/lib/solid_earth_tides)

# Add the costg directory conditionally
if(BUILD_COSTG)
  message(STATUS "Building costg binaries...")
  add_subdirectory(costg-benchmark/bin)
endif()


# The tests
if(BUILD_TESTING)
  include(CTest)
  add_subdirectory(test/unit_tests)
  add_subdirectory(test/should_not_compile)
  find_library(sofa sofa_c)
  if (sofa)
    add_subdirectory(test/sofa)
    message(STATUS "found sofa lib, will compile relevant tests")
  else()
    message(STATUS "sofa lib not found; tests in test/sofa will not be compiled.")
  endif()
  enable_testing()
endif()

# install library
install(TARGETS iers
        EXPORT iersTargets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
)

# Install public headers at: $PREFIX/iers/...
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/
	DESTINATION include/iers
)

# Install private headers to a separate directory $PREFIX/iers/core
install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/core/ 
        DESTINATION include/iers/core
)

install(EXPORT iersTargets
        FILE iersTargets.cmake
        NAMESPACE dso::
        DESTINATION lib/cmake
)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  "iersConfigVersion.cmake"
  VERSION ${iers_version}
  COMPATIBILITY AnyNewerVersion
)

install(FILES "iersConfig.cmake" "${CMAKE_CURRENT_BINARY_DIR}/iersConfigVersion.cmake"
  DESTINATION lib/cmake
)
