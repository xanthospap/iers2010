target_sources(iers
  PRIVATE
    ${CMAKE_SOURCE_DIR}/src/lib/doodson/doodson_from_chars.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/doodson/doodson_to_str.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/doodson/match_tadal_constituents.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/doodson/resolve_iers2010_doodson_string.cpp
    ${CMAKE_SOURCE_DIR}/src/lib/doodson/tidal_wave_tables.cpp
)
