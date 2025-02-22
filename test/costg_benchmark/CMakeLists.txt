add_compile_options(-Wno-unused-but-set-variable)
add_compile_options(-Wno-unused-variable)

include_directories(${CMAKE_SOURCE_DIR}/costg-benchmark/bin)

# JPL Ephemeris file(s)
find_file(DE
  NAMES "de421.bsp"
  PATHS ${CMAKE_SOURCE_DIR}/data/test_data
  NO_DEFAULT_PATH
  REQUIRED
)
find_file(DE_TLS
  NAMES "latest_leapseconds.tls"
  PATHS ${CMAKE_SOURCE_DIR}/data/test_data
  NO_DEFAULT_PATH
  REQUIRED
)

add_executable(tcheck_02gravityfield_icrf check_02gravityfield_icrf.cpp ${CMAKE_SOURCE_DIR}/costg-benchmark/bin/costg_utils.cpp)
target_link_libraries(tcheck_02gravityfield_icrf PRIVATE iers ${PROJECT_DEPENDENCIES})
target_include_directories(tcheck_02gravityfield_icrf PRIVATE ${CMAKE_SOURCE_DIR}/costg-benchmark/bin)
add_test(NAME costg_02_icrf COMMAND tcheck_02gravityfield_icrf 
  ${COSTG_DATA_DIR}/satellite/00orbit_itrf.txt 
  ${COSTG_DATA_DIR}/satellite/01earthRotation_rotaryMatrix.txt 
  ${COSTG_DATA_DIR}/models/EIGEN6-C4.gfc 
  ${COSTG_DATA_DIR}/satellite/02gravityfield_icrf.txt
)

add_executable(tcheck_02gravityfield_itrf check_02gravityfield_itrf.cpp ${CMAKE_SOURCE_DIR}/costg-benchmark/bin/costg_utils.cpp)
target_link_libraries(tcheck_02gravityfield_itrf PRIVATE iers ${PROJECT_DEPENDENCIES})
target_include_directories(tcheck_02gravityfield_itrf PRIVATE ${CMAKE_SOURCE_DIR}/costg-benchmark/bin)
add_test(NAME costg_02_itrf COMMAND tcheck_02gravityfield_itrf 
  ${COSTG_DATA_DIR}/satellite/00orbit_itrf.txt 
  ${COSTG_DATA_DIR}/models/EIGEN6-C4.gfc 
  ${COSTG_DATA_DIR}/satellite/02gravityfield_itrf.txt
)

add_executable(tcheck_03directTideSun_icrf check_03directTideSun_icrf.cpp ${CMAKE_SOURCE_DIR}/costg-benchmark/bin/costg_utils.cpp)
target_link_libraries(tcheck_03directTideSun_icrf PRIVATE iers ${PROJECT_DEPENDENCIES})
target_include_directories(tcheck_03directTideSun_icrf PRIVATE ${CMAKE_SOURCE_DIR}/costg-benchmark/bin)
add_test(NAME costg_03_sun COMMAND tcheck_03directTideSun_icrf 
  ${COSTG_DATA_DIR}/satellite/00orbit_icrf.txt 
  ${COSTG_DATA_DIR}/satellite/03directTideSun_icrf.txt
  ${DE} ${DE_TLS}
)

add_executable(tcheck_03directTideMoon_icrf check_03directTideMoon_icrf.cpp ${CMAKE_SOURCE_DIR}/costg-benchmark/bin/costg_utils.cpp)
target_link_libraries(tcheck_03directTideMoon_icrf PRIVATE iers ${PROJECT_DEPENDENCIES})
target_include_directories(tcheck_03directTideMoon_icrf PRIVATE ${CMAKE_SOURCE_DIR}/costg-benchmark/bin)
add_test(NAME costg_03_moon COMMAND tcheck_03directTideMoon_icrf 
  ${COSTG_DATA_DIR}/satellite/00orbit_icrf.txt 
  ${COSTG_DATA_DIR}/satellite/03directTideMoon_icrf.txt
  ${DE} ${DE_TLS}
)

##
##  Warning !!
##  Eop file is taken from the (local) test_data folder, NOT the COST-G
##  distributed one.
##
add_executable(tcheck_04solidEarthTide_icrf check_04solidEarthTide_icrf.cpp ${CMAKE_SOURCE_DIR}/costg-benchmark/bin/costg_utils.cpp)
target_link_libraries(tcheck_04solidEarthTide_icrf PRIVATE iers ${PROJECT_DEPENDENCIES})
target_include_directories(tcheck_04solidEarthTide_icrf PRIVATE ${CMAKE_SOURCE_DIR}/costg-benchmark/bin)
add_test(NAME costg_04 COMMAND tcheck_04solidEarthTide_icrf 
  ${COSTG_DATA_DIR}/satellite/00orbit_itrf.txt 
  ${COSTG_DATA_DIR}/satellite/04solidEarthTide_icrf.txt
  ${DE} ${DE_TLS}
  ${CMAKE_SOURCE_DIR}/data/test_data/eopc04.1962-now
  ${COSTG_DATA_DIR}/satellite/01earthRotation_rotaryMatrix.txt
)

add_executable(tcheck_05poleTide_icrf check_05poleTide_icrf.cpp ${CMAKE_SOURCE_DIR}/costg-benchmark/bin/costg_utils.cpp)
target_link_libraries(tcheck_05poleTide_icrf PRIVATE iers ${PROJECT_DEPENDENCIES})
target_include_directories(tcheck_05poleTide_icrf PRIVATE ${CMAKE_SOURCE_DIR}/costg-benchmark/bin)
add_test(NAME costg_05 COMMAND tcheck_05poleTide_icrf 
  ${COSTG_DATA_DIR}/satellite/00orbit_itrf.txt 
  ${COSTG_DATA_DIR}/satellite/01earthRotation_rotaryMatrix.txt
  ${COSTG_DATA_DIR}/satellite/01earthRotation_interpolatedEOP.txt
  ${COSTG_DATA_DIR}/satellite/05poleTide_icrf.txt
)

add_executable(tcheck_06oceanPoleTide_icrf check_06oceanPoleTide_icrf.cpp ${CMAKE_SOURCE_DIR}/costg-benchmark/bin/costg_utils.cpp)
target_link_libraries(tcheck_06oceanPoleTide_icrf PRIVATE iers ${PROJECT_DEPENDENCIES})
target_include_directories(tcheck_06oceanPoleTide_icrf PRIVATE ${CMAKE_SOURCE_DIR}/costg-benchmark/bin)
add_test(NAME costg_06 COMMAND tcheck_06oceanPoleTide_icrf 
  ${COSTG_DATA_DIR}/satellite/00orbit_itrf.txt 
  ${COSTG_DATA_DIR}/satellite/01earthRotation_rotaryMatrix.txt
  ${COSTG_DATA_DIR}/satellite/01earthRotation_interpolatedEOP.txt
  ${COSTG_DATA_DIR}/satellite/06oceanPoleTide_icrf.txt
  ${COSTG_DATA_DIR}/models/desaiscopolecoef.txt
)

add_executable(tcheck_07relativistic_icrf check_07relativistic_icrf.cpp ${CMAKE_SOURCE_DIR}/costg-benchmark/bin/costg_utils.cpp)
target_link_libraries(tcheck_07relativistic_icrf PRIVATE iers ${PROJECT_DEPENDENCIES})
target_include_directories(tcheck_07relativistic_icrf PRIVATE ${CMAKE_SOURCE_DIR}/costg-benchmark/bin)
add_test(NAME costg_07 COMMAND tcheck_07relativistic_icrf 
  ${COSTG_DATA_DIR}/satellite/00orbit_icrf.txt 
  ${COSTG_DATA_DIR}/satellite/07relativistic_icrf.txt
  ${DE} ${DE_TLS}
)

add_executable(tcheck_08aod1b_RL06_icrf check_08aod1b_RL06_icrf.cpp ${CMAKE_SOURCE_DIR}/costg-benchmark/bin/costg_utils.cpp)
target_link_libraries(tcheck_08aod1b_RL06_icrf PRIVATE iers ${PROJECT_DEPENDENCIES})
target_include_directories(tcheck_08aod1b_RL06_icrf PRIVATE ${CMAKE_SOURCE_DIR}/costg-benchmark/bin)
add_test(NAME costg_08 COMMAND tcheck_08aod1b_RL06_icrf 
  ${COSTG_DATA_DIR}/satellite/00orbit_itrf.txt 
  ${COSTG_DATA_DIR}/models/AOD1B_2008-07-03_X_06.asc
  ${COSTG_DATA_DIR}/satellite/08aod1b_RL06_icrf.txt
  ${COSTG_DATA_DIR}/satellite/01earthRotation_rotaryMatrix.txt
  ${COSTG_DATA_DIR}/models
)

add_executable(tcheck_11oceanTide_fes2014b_34major_icrf check_11oceanTide_fes2014b_34major_icrf.cpp ${CMAKE_SOURCE_DIR}/costg-benchmark/bin/costg_utils.cpp)
target_link_libraries(tcheck_11oceanTide_fes2014b_34major_icrf PRIVATE iers ${PROJECT_DEPENDENCIES})
target_include_directories(tcheck_11oceanTide_fes2014b_34major_icrf PRIVATE ${CMAKE_SOURCE_DIR}/costg-benchmark/bin)
add_test(NAME costg_11_fes14b_major COMMAND tcheck_11oceanTide_fes2014b_34major_icrf 
  ${COSTG_DATA_DIR}/satellite/00orbit_itrf.txt 
  ${COSTG_DATA_DIR}/satellite/11oceanTide_fes2014b_34major_icrf.txt
  ${COSTG_DATA_DIR}/satellite/01earthRotation_rotaryMatrix.txt
  ${CMAKE_SOURCE_DIR}/data/test_data/eopc04.1962-now
  #${COSTG_DATA_DIR}/../oceanAndAtmosphericTides/models/FES2014b_OCN/FES2014b_OCN_001fileList.txt
  #${COSTG_DATA_DIR}/../oceanAndAtmosphericTides/models/FES2014b_OCN/
  ${COSTG_DATA_DIR}/models/FES2014b_oceanTide/fes2014b/FES2014b_001fileList.txt
  ${COSTG_DATA_DIR}/models/FES2014b_oceanTide/fes2014b
)

add_executable(tcheck_11oceanTide_fes2014b_with361Admittance_icrf check_11oceanTide_fes2014b_with361Admittance_icrf.cpp ${CMAKE_SOURCE_DIR}/costg-benchmark/bin/costg_utils.cpp)
target_link_libraries(tcheck_11oceanTide_fes2014b_with361Admittance_icrf PRIVATE iers ${PROJECT_DEPENDENCIES})
target_include_directories(tcheck_11oceanTide_fes2014b_with361Admittance_icrf PRIVATE ${CMAKE_SOURCE_DIR}/costg-benchmark/bin)
add_test(NAME costg_11_fes14b_admittance COMMAND tcheck_11oceanTide_fes2014b_with361Admittance_icrf 
  ${COSTG_DATA_DIR}/satellite/00orbit_itrf.txt 
  ${COSTG_DATA_DIR}/satellite/11oceanTide_fes2014b_with361Admittance_icrf.txt
  ${COSTG_DATA_DIR}/satellite/01earthRotation_rotaryMatrix.txt
  ${CMAKE_SOURCE_DIR}/data/test_data/eopc04.1962-now
  #${COSTG_DATA_DIR}/../oceanAndAtmosphericTides/models/FES2014b_OCN/FES2014b_OCN_001fileList.txt
  #${COSTG_DATA_DIR}/../oceanAndAtmosphericTides/models/FES2014b_OCN/
  ${COSTG_DATA_DIR}/models/FES2014b_oceanTide/fes2014b/FES2014b_001fileList.txt
  ${COSTG_DATA_DIR}/../oceanAndAtmosphericTides/models/FES2014b_OCN/FES2014b_OCN_002doodson.txt
  #${COSTG_DATA_DIR}/../oceanAndAtmosphericTides/models/FES2014b_OCN/FES2014b_OCN_003admittance.txt
  ${COSTG_DATA_DIR}/models/FES2014b_oceanTide/admittance/fes2014b_admittance_linear_linear.txt
  ${COSTG_DATA_DIR}/models/FES2014b_oceanTide/fes2014b
)

add_executable(tcheck_10oceanTide_eot11a_M2_icrf check_10oceanTide_eot11a_M2_icrf.cpp ${CMAKE_SOURCE_DIR}/costg-benchmark/bin/costg_utils.cpp)
target_link_libraries(tcheck_10oceanTide_eot11a_M2_icrf PRIVATE iers ${PROJECT_DEPENDENCIES})
target_include_directories(tcheck_10oceanTide_eot11a_M2_icrf PRIVATE ${CMAKE_SOURCE_DIR}/costg-benchmark/bin)
add_test(NAME costg_10_eot11a_m2 COMMAND tcheck_10oceanTide_eot11a_M2_icrf 
  ${COSTG_DATA_DIR}/satellite/00orbit_itrf.txt 
  ${COSTG_DATA_DIR}/satellite/10oceanTide_eot11a_M2_icrf.txt
  ${COSTG_DATA_DIR}/satellite/01earthRotation_rotaryMatrix.txt
  ${CMAKE_SOURCE_DIR}/data/test_data/eopc04.1962-now
  ${COSTG_DATA_DIR}/models/EOT11a_oceanTide/eot11a/eot11aM2_001fileList.txt
  ${COSTG_DATA_DIR}/models/EOT11a_oceanTide/eot11a
)

add_executable(tcheck_10oceanTide_eot11a_18major_icrf check_10oceanTide_eot11a_M2_icrf.cpp ${CMAKE_SOURCE_DIR}/costg-benchmark/bin/costg_utils.cpp)
target_link_libraries(tcheck_10oceanTide_eot11a_18major_icrf PRIVATE iers ${PROJECT_DEPENDENCIES})
target_include_directories(tcheck_10oceanTide_eot11a_18major_icrf PRIVATE ${CMAKE_SOURCE_DIR}/costg-benchmark/bin)
add_test(NAME costg_10_eot11a_major COMMAND tcheck_10oceanTide_eot11a_18major_icrf 
  ${COSTG_DATA_DIR}/satellite/00orbit_itrf.txt 
  ${COSTG_DATA_DIR}/satellite/10oceanTide_eot11a_18major_icrf.txt
  ${COSTG_DATA_DIR}/satellite/01earthRotation_rotaryMatrix.txt
  ${CMAKE_SOURCE_DIR}/data/test_data/eopc04.1962-now
  ${COSTG_DATA_DIR}/models/EOT11a_oceanTide/eot11a/eot11a_001fileList.txt
  ${COSTG_DATA_DIR}/models/EOT11a_oceanTide/eot11a
)

add_executable(tcheck_10oceanTide_eot11a_with256Admittance_icrf check_10oceanTide_eot11a_with256Admittance_icrf.cpp ${CMAKE_SOURCE_DIR}/costg-benchmark/bin/costg_utils.cpp)
target_link_libraries(tcheck_10oceanTide_eot11a_with256Admittance_icrf PRIVATE iers ${PROJECT_DEPENDENCIES})
target_include_directories(tcheck_10oceanTide_eot11a_with256Admittance_icrf PRIVATE ${CMAKE_SOURCE_DIR}/costg-benchmark/bin)
add_test(NAME costg_10_eot11a_admittance COMMAND tcheck_10oceanTide_eot11a_with256Admittance_icrf 
  ${COSTG_DATA_DIR}/satellite/00orbit_itrf.txt 
  ${COSTG_DATA_DIR}/satellite/10oceanTide_eot11a_with256Admittance_icrf.txt
  ${COSTG_DATA_DIR}/satellite/01earthRotation_rotaryMatrix.txt
  ${CMAKE_SOURCE_DIR}/data/test_data/eopc04.1962-now
  ${COSTG_DATA_DIR}/models/EOT11a_oceanTide/eot11a/eot11a_001fileList.txt
  ${COSTG_DATA_DIR}/models/EOT11a_oceanTide/eot11a/EOT11A_002doodson.txt
  ${COSTG_DATA_DIR}/models/EOT11a_oceanTide/eot11a/EOT11A_003admittance.txt
  ${COSTG_DATA_DIR}/models/EOT11a_oceanTide/eot11a
)
