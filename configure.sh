#! /bin/bash

##
##  This scripts is a helper for installing the libiers10++ library.
##

## Check for gcc and its version
if ! gcc --version 2>/dev/null
then
    echo "Could not find the gnu gcc compiler!"
    exit 1
fi


