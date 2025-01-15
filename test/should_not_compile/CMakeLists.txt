add_compile_options(-Wno-unused-but-set-variable)
add_compile_options(-Wno-unused-variable)

# --------------------------------------------------------------------
# Test coeff_matrix_2d_rowise_mock
# --------------------------------------------------------------------
add_executable(coeff_matrix_2d_rowise_mock coeff_matrix_2d_rowise_mock.cpp)
target_link_libraries(coeff_matrix_2d_rowise_mock PRIVATE iers)
target_link_libraries(coeff_matrix_2d_rowise_mock PRIVATE ${PROJECT_DEPENDENCIES})
set_target_properties(coeff_matrix_2d_rowise_mock PROPERTIES
    EXCLUDE_FROM_ALL TRUE
    EXCLUDE_FROM_DEFAULT_BUILD TRUE
)
add_test(
    NAME coeff_matrix_2d_rowise_mock
    COMMAND ${CMAKE_COMMAND} --build . --target coeff_matrix_2d_rowise_mock --config $<CONFIG>
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(coeff_matrix_2d_rowise_mock PROPERTIES WILL_FAIL TRUE)

# --------------------------------------------------------------------
# Test coeff_matrix_2d_tricolwise22_mock
# --------------------------------------------------------------------
add_executable(coeff_matrix_2d_tricolwise22_mock coeff_matrix_2d_tricolwise22_mock.cpp)
target_link_libraries(coeff_matrix_2d_tricolwise22_mock PRIVATE iers)
target_link_libraries(coeff_matrix_2d_tricolwise22_mock PRIVATE ${PROJECT_DEPENDENCIES})
set_target_properties(coeff_matrix_2d_tricolwise22_mock PROPERTIES
    EXCLUDE_FROM_ALL TRUE
    EXCLUDE_FROM_DEFAULT_BUILD TRUE
)
add_test(
    NAME coeff_matrix_2d_tricolwise22_mock
    COMMAND ${CMAKE_COMMAND} --build . --target coeff_matrix_2d_tricolwise22_mock --config $<CONFIG>
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(coeff_matrix_2d_tricolwise22_mock PROPERTIES WILL_FAIL TRUE)

# --------------------------------------------------------------------
# Test coeff_matrix_2d_tricolwise_cp_mock
# --------------------------------------------------------------------
#add_executable(coeff_matrix_2d_tricolwise_cp_mock coeff_matrix_2d_tricolwise_cp_mock.cpp)
#target_link_libraries(coeff_matrix_2d_tricolwise_cp_mock PRIVATE iers)
#target_link_libraries(coeff_matrix_2d_tricolwise_cp_mock PRIVATE ${PROJECT_DEPENDENCIES})
#set_target_properties(coeff_matrix_2d_tricolwise_cp_mock PROPERTIES
#    EXCLUDE_FROM_ALL TRUE
#    EXCLUDE_FROM_DEFAULT_BUILD TRUE
#)
#add_test(
#    NAME coeff_matrix_2d_tricolwise_cp_mock
#    COMMAND ${CMAKE_COMMAND} --build . --target coeff_matrix_2d_tricolwise_cp_mock --config $<CONFIG>
#    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#)
#set_tests_properties(coeff_matrix_2d_tricolwise_cp_mock PROPERTIES WILL_FAIL TRUE)

# --------------------------------------------------------------------
# Test coeff_matrix_2d_trirowise_mock
# --------------------------------------------------------------------
add_executable(coeff_matrix_2d_trirowise_mock coeff_matrix_2d_trirowise_mock.cpp)
target_link_libraries(coeff_matrix_2d_trirowise_mock PRIVATE iers)
target_link_libraries(coeff_matrix_2d_trirowise_mock PRIVATE ${PROJECT_DEPENDENCIES})
set_target_properties(coeff_matrix_2d_trirowise_mock PROPERTIES
    EXCLUDE_FROM_ALL TRUE
    EXCLUDE_FROM_DEFAULT_BUILD TRUE
)
add_test(
    NAME coeff_matrix_2d_trirowise_mock
    COMMAND ${CMAKE_COMMAND} --build . --target coeff_matrix_2d_trirowise_mock --config $<CONFIG>
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(coeff_matrix_2d_trirowise_mock PROPERTIES WILL_FAIL TRUE)

# --------------------------------------------------------------------
# Test coeff_matrix_2d_colwise_mock
# --------------------------------------------------------------------
add_executable(coeff_matrix_2d_colwise_mock coeff_matrix_2d_colwise_mock.cpp)
target_link_libraries(coeff_matrix_2d_colwise_mock PRIVATE iers)
target_link_libraries(coeff_matrix_2d_colwise_mock PRIVATE ${PROJECT_DEPENDENCIES})
set_target_properties(coeff_matrix_2d_colwise_mock PROPERTIES
    EXCLUDE_FROM_ALL TRUE
    EXCLUDE_FROM_DEFAULT_BUILD TRUE
)
add_test(
    NAME coeff_matrix_2d_colwise_mock
    COMMAND ${CMAKE_COMMAND} --build . --target coeff_matrix_2d_colwise_mock --config $<CONFIG>
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(coeff_matrix_2d_colwise_mock PROPERTIES WILL_FAIL TRUE)

# --------------------------------------------------------------------
# Test coeff_matrix_2d_tricolwise21_mock
# --------------------------------------------------------------------
add_executable(coeff_matrix_2d_tricolwise21_mock coeff_matrix_2d_tricolwise21_mock.cpp)
target_link_libraries(coeff_matrix_2d_tricolwise21_mock PRIVATE iers)
target_link_libraries(coeff_matrix_2d_tricolwise21_mock PRIVATE ${PROJECT_DEPENDENCIES})
set_target_properties(coeff_matrix_2d_tricolwise21_mock PROPERTIES
    EXCLUDE_FROM_ALL TRUE
    EXCLUDE_FROM_DEFAULT_BUILD TRUE
)
add_test(
    NAME coeff_matrix_2d_tricolwise21_mock
    COMMAND ${CMAKE_COMMAND} --build . --target coeff_matrix_2d_tricolwise21_mock --config $<CONFIG>
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(coeff_matrix_2d_tricolwise21_mock PROPERTIES WILL_FAIL TRUE)

# --------------------------------------------------------------------
# Test coeff_matrix_2d_tricolwise_mock
# --------------------------------------------------------------------
add_executable(coeff_matrix_2d_tricolwise_mock coeff_matrix_2d_tricolwise_mock.cpp)
target_link_libraries(coeff_matrix_2d_tricolwise_mock PRIVATE iers)
target_link_libraries(coeff_matrix_2d_tricolwise_mock PRIVATE ${PROJECT_DEPENDENCIES})
set_target_properties(coeff_matrix_2d_tricolwise_mock PROPERTIES
    EXCLUDE_FROM_ALL TRUE
    EXCLUDE_FROM_DEFAULT_BUILD TRUE
)
add_test(
    NAME coeff_matrix_2d_tricolwise_mock
    COMMAND ${CMAKE_COMMAND} --build . --target coeff_matrix_2d_tricolwise_mock --config $<CONFIG>
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_tests_properties(coeff_matrix_2d_tricolwise_mock PROPERTIES WILL_FAIL TRUE)
