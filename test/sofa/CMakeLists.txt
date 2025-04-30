add_compile_options(-Wno-unused-but-set-variable)
add_compile_options(-Wno-unused-variable)

# Define a macro to add the private include directory for all targets
macro(add_internal_includes target_name)
  target_include_directories(${target_name} PRIVATE ${CMAKE_SOURCE_DIR}/src)
endmacro()

set(TEST_SOURCES
  rot_pom.cpp
  rot_pn.cpp
  era00.cpp
  fundarg.cpp
  s06.cpp
  xy06a.cpp
  gmst.cpp
  gmst82.cpp
  sp00.cpp
  xys06fargs.cpp
)

# Process each source file and create an executable
foreach(SOURCE_FILE IN LISTS TEST_SOURCES)
  # Get the filename without an extension
  get_filename_component(EXECUTABLE_NAME ${SOURCE_FILE} NAME_WE)

  # Define the executable
  add_executable(${EXECUTABLE_NAME} ${SOURCE_FILE})

  # Link the executable to the required libraries
  target_link_libraries(${EXECUTABLE_NAME} PRIVATE iers ${PROJECT_DEPENDENCIES} ${SOFAC})

  target_include_directories(${EXECUTABLE_NAME}
    PRIVATE ${CMAKE_SOURCE_DIR}/src)

  add_test(NAME ${EXECUTABLE_NAME} COMMAND ${EXECUTABLE_NAME})
endforeach()
