#! /bin/bash

if ! test -f configure.ac
then
    echo "ERROR! Failed to find configure.ac"
    exit 1
fi    

echo "Running automake to create missing scripts"
automake --add-missing --copy

./configure