#include <stdio.h>
#include "iers2010.hpp"

int main (int argc,char *argv[])
{
    int status = 0;

    if (argc!=2) {
        printf ("\nUsage: test-hardisp <blq input file>\n");
        return 1;
    }

    printf ("\nTesting Function HARDISP");

    int idate[] = {2009, 6, 25, 0, 0, 0};
    int it_size = 6;
    int irnt = 24;
    double samp = 3600.0e0;
    
    status = iers2010::hardisp (idate,it_size,irnt,samp,argv[1]);
    
    if (status)
        printf ("\nError running function iers2010::hardisp");
    
    printf ("\n");
    return status;
}
