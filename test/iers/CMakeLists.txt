add_compile_options(-Wno-unused-but-set-variable)
add_compile_options(-Wno-unused-variable)

add_executable(test_ortho_eop eop_oceantide.cpp)
target_link_libraries(test_ortho_eop PRIVATE iers ${PROJECT_DEPENDENCIES})
add_test(NAME test_ortho_eop COMMAND test_ortho_eop 
  ${CMAKE_CURRENT_SOURCE_DIR}/data/ortho_eop.dat
)