add_compile_options(-Wno-unused-but-set-variable)
add_compile_options(-Wno-unused-variable)

# Define a macro to add the private include directory for all targets
macro(add_internal_includes target_name)
  target_include_directories(${target_name} PRIVATE ${CMAKE_SOURCE_DIR}/src)
endmacro()

set(NO_ARGS_SOURCES 
  coeff_matrix_2d_colwise_cp.cpp
  coeff_matrix_2d_colwise.cpp
  coeff_matrix2d_colwise_reduce_2.cpp
  coeff_matrix2d_colwise_reduce_3.cpp
  coeff_matrix2d_colwise_reduce.cpp
  coeff_matrix_2d_rowise.cpp
  coeff_matrix_2d_rowwise_cp.cpp
  coeff_matrix2d_rowwise_reduce_2.cpp
  coeff_matrix2d_rowwise_reduce.cpp
  coeff_matrix_2d_tricolwise2.cpp
  coeff_matrix_2d_tricolwise_cp.cpp
  coeff_matrix_2d_tricolwise.cpp
  coeff_matrix2d_tricolwise_reduce_2.cpp
  coeff_matrix2d_tricolwise_reduce.cpp
  coeff_matrix_2d_trirowise.cpp
  coeff_matrix_2d_trirowwise_cp.cpp
  coeff_matrix2d_trirowwise_reduce.cpp
  stokes2anm.cpp
  anm2stokes.cpp
  alfs1.cpp
  #itrs2gcrs_velocity.cpp
  c2i_vel.cpp
  make_dehanttideinel_input.cpp
  solid_earth_tide_displacement.cpp
)

# Process each source file and create an executable
foreach(SOURCE_FILE IN LISTS NO_ARGS_SOURCES)
  # Get the filename without an extension
  get_filename_component(EXECUTABLE_NAME ${SOURCE_FILE} NAME_WE)

  # Define the executable
  add_executable(${EXECUTABLE_NAME} ${SOURCE_FILE})

  # Link the executable to the required libraries
  target_link_libraries(${EXECUTABLE_NAME} PRIVATE iers ${PROJECT_DEPENDENCIES})

  target_include_directories(${EXECUTABLE_NAME}
    PRIVATE ${CMAKE_SOURCE_DIR}/src)

  add_test(NAME ${EXECUTABLE_NAME} COMMAND ${EXECUTABLE_NAME})
endforeach()

set(NO_ARGS_SOURCES_EXPFAIL
  coeff_matrix2d_tricolwise_reduce_runtime_error.cpp
  coeff_matrix2d_trirowwise_reduce_runtime_error.cpp
  coeff_matrix_2d_tricolwise_cp_runtime_error.cpp
  coeff_matrix_2d_tricolwise23_runtime_error.cpp
)

# Process each source file and create an executable
foreach(SOURCE_FILE IN LISTS NO_ARGS_SOURCES_EXPFAIL)
  # Get the filename without an extension
  get_filename_component(EXECUTABLE_NAME ${SOURCE_FILE} NAME_WE)

  # Define the executable
  add_executable(${EXECUTABLE_NAME} ${SOURCE_FILE})

  # Link the executable to the required libraries
  target_link_libraries(${EXECUTABLE_NAME} PRIVATE iers ${PROJECT_DEPENDENCIES})

  target_include_directories(${EXECUTABLE_NAME}
    PRIVATE ${CMAKE_SOURCE_DIR}/src)

  add_test(NAME ${EXECUTABLE_NAME} 
    COMMAND ${CMAKE_COMMAND} -E env $<TARGET_FILE:${EXECUTABLE_NAME}>)

  # Specify that the test is expected to fail
  set_tests_properties(${EXECUTABLE_NAME} PROPERTIES 
    WILL_FAIL TRUE
    FAIL_REGULAR_EXPRESSION "Aborted \\(core dumped\\)|rows\\(\\) >= cprows && cols\\(\\) >= cpcols|SIGABRT|Assertion .* failed"
  )
endforeach()

##
##  Testing EOPs :
##   
##  (1) & (2) Test parsing C04-14 and C04-20 series
##  (3)       Test interpolation; actual interpolation values are not checked, 
##            only responses, status, etc ...
##
add_executable(eop1 eop1.cpp)
target_link_libraries(eop1 PRIVATE iers ${PROJECT_DEPENDENCIES})
target_include_directories(eop1 PRIVATE ${CMAKE_SOURCE_DIR}/src)
find_file(EOPC0420
  NAMES "eopc04.1962-now"
  PATHS ${CMAKE_SOURCE_DIR}/data/test_data
  NO_DEFAULT_PATH
  REQUIRED)
add_test(NAME eopc0420 COMMAND eop1 ${EOPC0420})
find_file(EOPC0414
  NAMES "eopc04_IAU2000.62-now"
  PATHS ${CMAKE_SOURCE_DIR}/data/test_data
  NO_DEFAULT_PATH
  REQUIRED)
add_test(NAME eopc0414 COMMAND eop1 ${EOPC0420})
add_executable(eop2 eop2.cpp)
target_link_libraries(eop2 PRIVATE iers ${PROJECT_DEPENDENCIES})
target_include_directories(eop2 PRIVATE ${CMAKE_SOURCE_DIR}/src)
add_test(NAME eop_interpolation COMMAND eop2 ${EOPC0420})

##
## Tests with AOD1B non-tidal product file
##
add_executable(test_aod1bin_collect_coeffs_1 test_aod1bin_collect_coeffs_1.cpp)
target_link_libraries(test_aod1bin_collect_coeffs_1 PRIVATE iers ${PROJECT_DEPENDENCIES})
target_include_directories(test_aod1bin_collect_coeffs_1 PRIVATE ${CMAKE_SOURCE_DIR}/src)
find_file(AOD1B
  NAMES "AOD1B_2008-07-03_X_06.asc"
  PATHS ${CMAKE_SOURCE_DIR}/data/test_data
  NO_DEFAULT_PATH
  REQUIRED)
add_test(NAME aod1b_01 COMMAND test_aod1bin_collect_coeffs_1 ${AOD1B})

##
## Tests with ocean loading models
##
add_executable(test_ocean_tides_1 ocean_tides_instances.cpp)
target_link_libraries(test_ocean_tides_1 PRIVATE iers ${PROJECT_DEPENDENCIES})
target_include_directories(test_ocean_tides_1 PRIVATE ${CMAKE_SOURCE_DIR}/src)
find_path(OCEAN_TIDES_DIR
  NAMES test_ocean_tide_models.md
  PATHS ${CMAKE_SOURCE_DIR}/data/test_data/ocean_tides
  NO_DEFAULT_PATH
  REQUIRED)
add_test(NAME make_ocean_tides COMMAND test_ocean_tides_1 ${OCEAN_TIDES_DIR})

##
## Terrestrial-to-Celestial Vs Astropy
##  
find_file(EOPN
  NAMES eopc0420.b2000e2023
  PATHS ${CMAKE_SOURCE_DIR}/data/test_data
  NO_DEFAULT_PATH
  REQUIRED)

add_executable(i2c_astropy i2c_astropy.cpp)
target_link_libraries(i2c_astropy PRIVATE iers ${PROJECT_DEPENDENCIES})
target_include_directories(i2c_astropy PRIVATE ${CMAKE_SOURCE_DIR}/src)

## sat::swot
find_file(SWOT
  NAMES swot_astropy_ref
  PATHS ${CMAKE_SOURCE_DIR}/data/test_data
  NO_DEFAULT_PATH
  REQUIRED)
add_test(NAME i2c_astropy_swot COMMAND i2c_astropy ${SWOT} ${EOPN})
## sat::cryosat-2
find_file(CS2
  NAMES cs2_astropy_ref
  PATHS ${CMAKE_SOURCE_DIR}/data/test_data
  NO_DEFAULT_PATH
  REQUIRED)
add_test(NAME i2c_astropy_cs2 COMMAND i2c_astropy ${CS2} ${EOPN})
## sat::jason-3
find_file(JA3
  NAMES ja3_astropy_ref
  PATHS ${CMAKE_SOURCE_DIR}/data/test_data
  NO_DEFAULT_PATH
  REQUIRED)
add_test(NAME i2c_astropy_ja3 COMMAND i2c_astropy ${JA3} ${EOPN})
