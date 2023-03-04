#! /usr/bin/bash

## path to DIR
RPTH="${HOME}/Software/iers2010/test/sofa_unit_tests"

## list of programs to run (assuming .out extension)
progs=("sofa-fundamental-angles" "sofa-test-ee06a" "sofa-test-era00" "sofa-test-gmst00" "sofa-test-gmst06")

## make sure all programs exist
for p in ${progs[@]}; do
  if ! test -f ${RPTH}/${p}.out ; then
    echo "Failed to find program ${RPTH}/${p}.out"
    exit 1
  fi
done

err=0
## loop through and execute; store status
for p in ${progs[@]}; do
  (${RPTH}/${p}.out)
  err=$(($err + $?))
done

exit $err
