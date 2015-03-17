#! /bin/bash

awk 'NR==FNR {
        for (i=1;i<=NF;i++) 
            a[FNR,i]=$i;
        }
    NF {
                printf ("%+010.7f %+010.7f %+010.7f\n",a[FNR,1]-$1,a[FNR,2]-$2,a[FNR,3]-$3)
    }
    ' ${1} ${2}
