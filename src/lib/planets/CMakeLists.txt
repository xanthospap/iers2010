target_sources(iers
  PRIVATE
    ${CMAKE_SOURCE_DIR}/src/lib/planets/kernel_consts.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/planets/kernel.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/planets/planet_pos.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/planets/planet_state.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/planets/planet_to_naif_id.cpp
)
