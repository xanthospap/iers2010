#! /usr/bin/bash

set -e

## Source code. Links to IERS2010 website.
links=("https://iers-conventions.obspm.fr/content/chapter8/software/ORTHO_EOP.F" \
"https://iers-conventions.obspm.fr/content/chapter8/software/CNMTX.F" \
"https://iers-conventions.obspm.fr/content/chapter8/software/RG_ZONT2.F" \
"https://iers-conventions.obspm.fr/content/chapter8/software/FUNDARG.F" \
"https://iers-conventions.obspm.fr/content/chapter5/software/UTLIBR.F" \
"https://iers-conventions.obspm.fr/content/chapter5/software/PMSDNUT2.F")

## Download every link (if basename does not exist).
for url in "${links[@]}"; do 
  filename=$(basename "$url")
  if [[ ! -f "$filename" ]]; then
    wget "$url" 
  fi
done

## Run the makefile.
make
