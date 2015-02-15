#! /bin/bash

DIR=/home/xanthos/myrepos/iers2010

for i in autom4te.cache \
     aclocal.m4 \
     config.h* \
     configure \
     Makefile \
     src/Makefile
do
    rm -rf ${DIR}/${i} 2>/dev/null
done
