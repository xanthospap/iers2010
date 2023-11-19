
data_files = [
    {"url":"http://icgem.gfz-potsdam.de/getmodel/gfc/4cc119d62d3f83adce857914bcabdfa7ca2f91677ed2905934cf52698584b4c9/EIGEN-6S4%20(v2).gfc", 
    "local": "data/EIGEN-6S4(v2).gfc"},
    {"url":"http://icgem.gfz-potsdam.de/getmodel/gfc/96e984f23a31505411d943341767e9ff082ad2c77fbbdcba1be2e393cf35110a/EIGEN-GRGS.RL04.MEAN-FIELD.gfc",
    "local": "data/EIGEN-GRGS.RL04.MEAN-FIELD.gfc"}
]

str1 = """Icgem filename: data/EIGEN-6S4v2.gfc
Product Type  : gravity_field
Model Name    : EIGEN-6S4v2
Tide System   : tide_free
Errors        : calibrated (sigma calibration factor =  2.00)
GM            : 398600441500000.000[m^3s^-2]
R_earth       : 6378136.459999999962747[m]
Max degree    : 300
"""

str2 = """Icgem filename: data/EIGEN-GRGS.RL04.MEAN-FIELD.gfc
Product Type  : gravity_field
Model Name    : EIGEN-GRGS.RL04.MEAN-FIELD.linear_mean_pole
Tide System   : tide_free
Errors        : formal (sigma calibration factor =  1.00)
GM            : 398600441500000.000[m^3s^-2]
R_earth       : 6378136.459999999962747[m]
Max degree    : 300
"""

special_progs = [
    {"prog":"parse-icgem-header.out", 
     "path": "test/unit_tests",
    "args": ["data/EIGEN-GRGS.RL04.MEAN-FIELD.gfc"], 
    "results": str2, "exit": 0},
    {"prog":"parse-icgem-header.out", 
     "path": "test/unit_tests",
    "args": ["data/EIGEN-6S4v2.gfc"], 
    "results": str1, "exit": 0}
]

