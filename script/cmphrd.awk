#! /bin/awk -f

##
##  This awk script, will compare two files, both holding 3 columns and
##+ the same number of lines and output their differenes. It is meant to
##+ be used for comparing output records from the hardisp routine.
##

function fabs(num)
{
    return (num>0) ? (num) : (num*-1);
}

BEGIN {
    printf ("------------------------------------------------------\n");
    printf (" # Line        Delta du      Delta dw        Delta ds (meters)\n")
    printf ("------------------------------------------------------\n");
    mean_u = 0.0e0;
    mean_w = 0.0e0;
    mean_s = 0.0e0;
}
NR==FNR {
    for (i=1;i<=NF;i++)
        a[FNR][i]=$i
    next
}
NF==3 {
    printf ("[LINE=%2i] %14.6f %14.6f %14.6f\n",FNR,
            fabs(a[FNR][1]-$1),fabs(a[FNR][2]-$2),fabs(a[FNR][3]-$3))
    mean_u += (a[FNR][1]-$1)
    mean_w += (a[FNR][2]-$2)
    mean_s += (a[FNR][3]-$3)
}
END {
    mean_u = mean_u / FNR
    printf ("------------------------------------------------------\n")
    printf ("AVERAGE : %14.6f %14.6f %14.6f\n",mean_u,mean_w,mean_s);
    printf ("------------------------------------------------------\n")
    printf ("          du : Radial tidal ocean loading displacement\n")
    printf ("          dw : West tidal ocean loading displacement\n")
    printf ("          ds : South tidal ocean loading displacement\n")
    printf ("------------------------------------------------------\n")
}
