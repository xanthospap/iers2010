#! /bin/bash

set -e

for c in c2t c2ixys era00 nut00a nut06a pfw06 rotation_matrix_3 s06 xy06
do
    g++ -Wall -std=c++17 -Wextra -Werror -pedantic -Wshadow -Winline -march=native -c ${c}.cpp -o ${c}.o
done

c=test_t2c
g++ -Wall -std=c++17 -Wextra -Werror -pedantic -Wshadow -Winline -march=native *.o ${c}.cc
