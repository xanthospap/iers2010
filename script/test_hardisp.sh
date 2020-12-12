#! /bin/bash

echoerr() { echo "$@" 1>&2; }

C_HARDISP="../src/hardisp"
F_HARDISP="../fortran_impl/HARDISP/hardisp"
TEST_BLQ=("../test/onsa.blq" "../test/reyk.blq")
DATE_STR="2009 6 25 1 10 45"
SAMPLE_NR=(1 10 24 124 347 586 600 601 891 1199 1200 1201 1960 1799 1800 1801)
SAMPLE_NT=(1800) 

## check that the programs are there and compiled
if ! test -f ${C_HARDISP} ; then
  echoerr "[ERROR] No C++ executable found: ${C_HARDISP}"
  exit 1
fi

if ! test -f ${F_HARDISP} ; then
  echoerr "[ERROR] No FORTRAN executable found: ${F_HARDISP}"
  exit 1
fi

for tf in ${TEST_BLQ[@]} ; do if ! test -f $tf ; then 
  echoerr "[ERROR] No test BLQ file found: $tf"
  exit 1
fi ; done

for sn in ${SAMPLE_NR[@]} ; do
  for si in ${SAMPLE_NT[@]} ; do
    for sf in ${TEST_BLQ[@]} ; do
      ${C_HARDISP} ${DATE_STR} ${sn} ${si} < ${sf} | sed '/^$/d' > test.c || { echo '${C_HARDISP} ${DATE_STR} ${sn} ${si} failed' ; exit 1; }
      #echo "${C_HARDISP} ${DATE_STR} ${sn} ${si} < ${sf} > test.c"
      ${F_HARDISP} ${DATE_STR} ${sn} ${si} < ${sf} > test.f || { echo '${F_HARDISP} ${DATE_STR} ${sn} ${si} failed' ; exit 1; }
      #echo "${F_HARDISP} ${DATE_STR} ${sn} ${si} < ${sf} > test.f"
      #echo "Max diffs for ${sn} ${si} :"
      paste test.c test.f | awk '{printf "%9.6f %9.6f %9.6f\n", $1-$4, $2-$5, $3-$6}' > paste.out
      awk 'BEGIN{c1=0e0; c2=0e0; c3=0e0} \
      {if (sqrt($1*$1)>c1) {c1=$1} \
        if (sqrt($2*$2)>c2) {c2=$2} \
        if (sqrt($3*$3)>c3) {c3=$3} } \
        END {printf "%10.6f %10.6f %10.6f\n", c1, c2, c3}' paste.out
    done
  done
done
