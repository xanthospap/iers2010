#! /usr/bin/bash

set -e

## Source code. Links to IERS2010 website.
links=("https://iers-conventions.obspm.fr/content/chapter8/software/ORTHO_EOP.F" "https://iers-conventions.obspm.fr/content/chapter8/software/CNMTX.F")

## Download every link.
for url in "${links[@]}"; do wget "$url"; done

## Run the makefile.
