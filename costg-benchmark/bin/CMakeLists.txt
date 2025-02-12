# costg-benchmark/bin/CMakeLists.txt

set(COSTG_SOURCES
  check_03directTideMoon_icrf.cpp
  check_03directTideSun_icrf.cpp
  check_04solidEarthTide_icrf.cpp
  check_05poleTide_icrf.cpp
  check_06oceanPoleTide_icrf.cpp
  check_07relativistic_icrf.cpp
  check_08aod1b_RL06_icrf.cpp
  check_09aod1b_atmosphericTides_icrf.cpp
  check_09aod1b_atmosphericTides_S1_icrf.cpp
#  check_11oceanTide_fes2014b_34major_icrf.cpp
)

# Define the shared utility source file
set(COSTG_UTILS_SOURCE costg_utils.cpp)

# Process each source file and create an executable
foreach(SOURCE_FILE IN LISTS COSTG_SOURCES)
  # Get the filename without an extension
  get_filename_component(EXECUTABLE_NAME ${SOURCE_FILE} NAME_WE)

  # Define the executable
  add_executable(${EXECUTABLE_NAME} ${SOURCE_FILE} ${COSTG_UTILS_SOURCE})

  # Include the costg directory for headers
  target_include_directories(${EXECUTABLE_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

  # Link the executable to the required libraries
  target_link_libraries(${EXECUTABLE_NAME} PRIVATE iers ${PROJECT_DEPENDENCIES})
endforeach()
