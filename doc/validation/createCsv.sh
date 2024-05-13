#! /usr/bin/env bash

set -e

IDIR=/home/xanthos/Software/iers2010
BDIR=${IDIR}/costg-benchmark/bin
CDIR=${IDIR}/costg-benchmark/softwareComparison
#CDIR=${IDIR}/costg-benchmark/outgoing/ITSG/COST-G/softwareComparison
ODIR=${IDIR}/doc/validation

${BDIR}/check-02gravityfield-itrf.out \
$CDIR/satellite/00orbit_itrf.txt \
$CDIR/models/EIGEN6-C4.gfc \
$CDIR/satellite/02gravityfield_itrf.txt > ${ODIR}/02gravityfield_itrf.csv

${BDIR}/check-02gravityfield-icrf.out \
$CDIR/satellite/00orbit_itrf.txt \
$CDIR/satellite/01earthRotation_rotaryMatrix.txt \
$CDIR/models/EIGEN6-C4.gfc \
$CDIR/satellite/02gravityfield_icrf.txt > ${ODIR}/02gravityfield_icrf.csv

${BDIR}/check-03directTideMoon-icrf.out \
$CDIR/satellite/00orbit_icrf.txt \
$CDIR/satellite/03directTideMoon_icrf.txt \
${IDIR}/data/de421.bsp ${IDIR}/data/latest_leapseconds.tls > ${ODIR}/03directTideMoon_icrf.csv

${BDIR}/check-03directTideSun-icrf.out \
$CDIR/satellite/00orbit_icrf.txt \
$CDIR/satellite/03directTideSun_icrf.txt \
${IDIR}/data/de421.bsp ${IDIR}/data/latest_leapseconds.tls > ${ODIR}/03directTideSun_icrf.csv

${BDIR}/check-04solidEarthTide-icrf.out \
$CDIR/satellite/00orbit_itrf.txt \
$CDIR/satellite/04solidEarthTide_icrf.txt \
${IDIR}/data/de421.bsp \
${IDIR}/data/latest_leapseconds.tls \
$CDIR/models/eopc04_14_IAU2000.62-now \
$CDIR/satellite/01earthRotation_rotaryMatrix.txt > ${ODIR}/04solidEarthTide_icrf.csv

${BDIR}/check-07relativistic-icrf.out \
$CDIR/satellite/00orbit_icrf.txt \
$CDIR/satellite/07relativistic_icrf.txt \
${IDIR}/data/de421.bsp \
${IDIR}/data/latest_leapseconds.tls > ${ODIR}/07relativistic_icrf.csv

${BDIR}/check-08aod1b-RL06-icrf.out \
$CDIR/satellite/00orbit_itrf.txt \
$CDIR/models/AOD1B_2008-07-03_X_06.asc \
$CDIR/satellite/08aod1b_RL06_icrf.txt \
$CDIR/satellite/01earthRotation_rotaryMatrix.txt \
$CDIR/models > ${ODIR}/08aod1b_RL06_icrf.csv

${BDIR}/check-09aod1b-atmosphericTides-S1-icrf.out \
$CDIR/satellite/00orbit_itrf.txt \
$CDIR/models/AOD1B_tides/AOD1B_ATM_S1_06.asc \
$CDIR/satellite/01earthRotation_rotaryMatrix.txt \
$CDIR/models/eopc04_14_IAU2000.62-now \
$CDIR/satellite/09aod1b_atmosphericTides_S1_icrf.txt > ${ODIR}/09aod1b_atmosphericTides_S1_icrf.csv
