#! /bin/bash

set -e

OPTL=2

CMP_FLAGS="-Wall -std=c++17 -Wextra -Werror -pedantic -Wshadow -Winline -march=native -O$OPTL"
for c in  c2ixys.cpp  c2t.cpp  era00.cpp  nut00a.cpp  nut06a.cpp  p06e.cpp  pfw06.cpp  \
 pn06.cpp  rotation_matrix_3.cpp  s00.cpp  s06.cpp  xy06.cpp eect00.cpp
do
    g++ $CMP_FLAGS -c ${c} -o ${c/.cpp/.o}
done

c=test_t2c
g++ $CMP_FLAGS *.o ${c}.cc -o ver1.out


CMP_FLAGS=${CMP_FLAGS}" -DUSE_CUSTOM_CTORS"
for c in  c2ixys.cpp  c2t.cpp  era00.cpp  nut00a.cpp  nut06a.cpp  p06e.cpp  pfw06.cpp  \
 pn06.cpp  rotation_matrix_3.cpp  s00.cpp  s06.cpp  xy06.cpp eect00.cpp
do
    g++ $CMP_FLAGS -c ${c} -o ${c/.cpp/.o}
done

c=test_t2c
g++ $CMP_FLAGS *.o ${c}.cc -o ver2.out