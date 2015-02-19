#! /bin/bash

for i in autom4te.cache \
    aclocal.m4 \
    config.* \
    configure \
    depcomp \
    install-sh \
    Makefile \
    Makefile.in \
    missing \
    stamp-h1 \
    test-driver
do
    rm -rf ../${i} ##2>/dev/null
done

for i in src test
do
    cd ../${i} ; rm Makefile Makefile.in 2>/dev/null ; cd ../script;
    cd ../${i} ; rm Makefile Makefile.in 2>/dev/null ; cd ../script;
done

rm ../libiers10*.a 2>/dev/null

rm -rf ../doc/Makefile ../doc/Makefile.in ../doc/api_html \
    ../doc/doxyfile.stamp ../doc/man ../Doxyfile.in 2>/dev/null
