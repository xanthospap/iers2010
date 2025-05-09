add_compile_options(-Wno-unused-but-set-variable)
add_compile_options(-Wno-unused-variable)

add_executable(test_ortho_eop eop_oceantide.cpp)
target_link_libraries(test_ortho_eop PRIVATE iers ${PROJECT_DEPENDENCIES})
add_test(NAME test_ortho_eop COMMAND test_ortho_eop 
  ${CMAKE_CURRENT_SOURCE_DIR}/data/ortho_eop.dat
)

add_executable(test_rgzont2 eop_rgzont2.cpp)
target_link_libraries(test_rgzont2 PRIVATE iers ${PROJECT_DEPENDENCIES})
add_test(NAME test_rgzont2 COMMAND test_rgzont2 
  ${CMAKE_CURRENT_SOURCE_DIR}/data/rg_zont2.dat
)

add_executable(test_libration eop_libration.cpp)
target_link_libraries(test_libration PRIVATE iers ${PROJECT_DEPENDENCIES})
add_test(NAME test_libration COMMAND test_libration 
  ${CMAKE_CURRENT_SOURCE_DIR}/data/libration.dat
)
