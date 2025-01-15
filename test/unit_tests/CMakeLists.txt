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
  alfs1.cpp
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
