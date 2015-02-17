#! /bin/bash

if ! test -f configure.ac
then
    echo "ERROR! Failed to find configure.ac"
    exit 1
fi    

echo "Running automake to create missing scripts"
automake --add-missing --copy

./configure

for i in aclocal.m4 \
    autom4te.cache
do
    rm -rf ${i} 2>/dev/null
done
