#! /bin/bash

for f in http://maia.usno.navy.mil/conv2010/convupdt/chapter5/PMSDNUT2.F \
    http://maia.usno.navy.mil/conv2010/chapter5/UTLIBR.F \
    http://maia.usno.navy.mil/conv2010/chapter5/FUNDARG.F \
    http://maia.usno.navy.mil/conv2010/convupdt/chapter5/FCNNUT.F \
    http://maia.usno.navy.mil/conv2010/convupdt/chapter7/ARG2.F \
    http://maia.usno.navy.mil/conv2010/convupdt/chapter7/IERS_CMP_2015.F \
    ftp://maia.usno.navy.mil/conv2010/convupdt/chapter7/dehanttideinel/CAL2JD.F \
    ftp://maia.usno.navy.mil/conv2010/convupdt/chapter7/dehanttideinel/DAT.F \
    ftp://maia.usno.navy.mil/conv2010/convupdt/chapter7/dehanttideinel/DEHANTTIDEINEL.F \
    ftp://maia.usno.navy.mil/conv2010/convupdt/chapter7/dehanttideinel/NORM8.F \
    ftp://maia.usno.navy.mil/conv2010/convupdt/chapter7/dehanttideinel/SPROD.F \
    ftp://maia.usno.navy.mil/conv2010/convupdt/chapter7/dehanttideinel/ST1IDIU.F \
    ftp://maia.usno.navy.mil/conv2010/convupdt/chapter7/dehanttideinel/ST1ISEM.F \
    ftp://maia.usno.navy.mil/conv2010/convupdt/chapter7/dehanttideinel/ST1L1.F \
    ftp://maia.usno.navy.mil/conv2010/convupdt/chapter7/dehanttideinel/STEP2DIU.F \
    ftp://maia.usno.navy.mil/conv2010/convupdt/chapter7/dehanttideinel/STEP2LON.F \
    ftp://maia.usno.navy.mil/conv2010/convupdt/chapter7/dehanttideinel/ZERO_VEC8.F \
    http://maia.usno.navy.mil/conv2010/chapter8/ORTHO_EOP.F \
    http://maia.usno.navy.mil/conv2010/chapter8/CNMTX.F \
    http://maia.usno.navy.mil/conv2010/convupdt/chapter8/RG_ZONT2.F
    do
        fr=$(basename $f)
        if ! test -f $fr
            then
                wget $f
        fi
done
