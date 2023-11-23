
data_files = [
    {"url":"http://icgem.gfz-potsdam.de/getmodel/gfc/4cc119d62d3f83adce857914bcabdfa7ca2f91677ed2905934cf52698584b4c9/EIGEN-6S4%20(v2).gfc", 
    "local": "data/EIGEN-6S4v2.gfc"},
    {"url":"http://icgem.gfz-potsdam.de/getmodel/gfc/96e984f23a31505411d943341767e9ff082ad2c77fbbdcba1be2e393cf35110a/EIGEN-GRGS.RL04.MEAN-FIELD.gfc",
    "local": "data/EIGEN-GRGS.RL04.MEAN-FIELD.gfc"},
    {"url": "http://icgem.gfz-potsdam.de/getmodel/gfc/fc51e0803271a9355bd55be4eea2790049c50a2e0a96a2077503da1df7db3bd0/EIGEN-5C.gfc",
    "local": "data/EIGEN-5C.gfc"}
]

str1 = """Icgem filename: data/EIGEN-6S4v2.gfc
Product Type  : gravity_field
Model Name    : EIGEN-6S4v2
Tide System   : tide_free
Errors        : calibrated
GM            : 398600441500000.000[m^3s^-2]
R_earth       : 6378136.459999999962747[m]
Max degree    : 300
"""

str2 = """Icgem filename: data/EIGEN-GRGS.RL04.MEAN-FIELD.gfc
Product Type  : gravity_field
Model Name    : EIGEN-GRGS.RL04.MEAN-FIELD.linear_mean_pole
Tide System   : tide_free
Errors        : formal
GM            : 398600441500000.000[m^3s^-2]
R_earth       : 6378136.459999999962747[m]
Max degree    : 300
"""

str3 = """Icgem filename: data/EIGEN-5C.gfc
Product Type  : gravity_field
Model Name    : EIGEN-5C
Tide System   : tide_free
Errors        : calibrated
GM            : 398600441500000.000[m^3s^-2]
R_earth       : 6378136.459999999962747[m]
Max degree    : 360
"""

str4 = """C = +9.048055032889907e-07 S = -6.189981297419080e-07
C = +9.048054845451114e-07 S = -6.189984387467481e-07
C = +9.048055422226470e-07 S = -6.189975126446038e-07
C = +9.048055225152388e-07 S = -6.189978210378097e-07
C = +9.048055032889907e-07 S = -6.189981297419080e-07
C = +9.048054845451114e-07 S = -6.189984387467481e-07
C = +9.048055422226470e-07 S = -6.189975126446038e-07
C = +9.048055225152388e-07 S = -6.189978210378097e-07
C = +9.047848403910832e-07 S = -6.189574273460872e-07
C = +9.047973678276128e-07 S = -6.189976324873710e-07
C = +9.047729241066045e-07 S = -6.189861582824630e-07
C = +9.048190638560238e-07 S = -6.190054537003117e-07
C = +9.047971428653616e-07 S = -6.190049132113741e-07
C = +9.047835455751224e-07 S = -6.190053202803755e-07
C = +9.047629841648657e-07 S = -6.190393959810343e-07
C = +9.047769826341694e-07 S = -6.190183462297375e-07
C = +9.047292131124359e-07 S = -6.190494668374820e-07
C = +9.047402184377304e-07 S = -6.189882384464848e-07
C = +9.047627311668372e-07 S = -6.190417353013939e-07
C = +9.047342376946078e-07 S = -6.190127540098040e-07
C = +9.047446052423377e-07 S = -6.190508232612063e-07
C = +9.047421030186882e-07 S = -6.190661746249831e-07
C = +9.047283335697546e-07 S = -6.190677084230348e-07
C = +9.047313873690279e-07 S = -6.190865078687803e-07
C = +9.047344419135161e-07 S = -6.191053076294447e-07
C = +9.047374972015834e-07 S = -6.191241077023384e-07
C = +9.047404375374159e-07 S = -6.191434697420685e-07
C = +9.047434913366892e-07 S = -6.191622691878140e-07
C = +9.047465458811774e-07 S = -6.191810689484784e-07
C = +9.047496011692447e-07 S = -6.191998690213720e-07
"""

special_progs = [
    {"prog":"parse-icgem-header.out", 
     "path": "test/unit_tests",
    "args": ["data/EIGEN-GRGS.RL04.MEAN-FIELD.gfc"], 
    "results": str2, "exit": 0},
    {"prog":"parse-icgem-header.out", 
     "path": "test/unit_tests",
    "args": ["data/EIGEN-6S4v2.gfc"], 
    "results": str1, "exit": 0},
    {"prog":"parse-icgem-header.out", 
     "path": "test/unit_tests",
    "args": ["data/EIGEN-5C.gfc"], 
    "results": str3, "exit": 0},
    {"prog":"parse-icgem-data.out", 
     "path": "test/unit_tests",
    "args": ["data/EIGEN-GRGS.RL04.MEAN-FIELD.gfc"], 
    "results": str4, "exit": 0}
]

