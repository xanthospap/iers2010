#include "iers2010.hpp"
#include "hardisp.hpp"
#include <stdexcept>
#include <iostream>
#include <fstream>

constexpr int ntin ( 11 );
double tamp[3][ntin], tph[3][ntin];

int collect_argv (std::ifstream& is)
{
    char line[256];
    char c = is.peek ();
    while (c == '$')
        is.getline (line,256);

    for (int i=0;i<3;i++) {
        for (int kk=0;kk<ntin;kk++) {
            is >> tamp[i][kk];
        }
        c = is.peek ();
        while (c == '$')
            is.getline (line,256);
    }
    if ( !is || is.fail() ) {
        printf ("\nError reading data lines. Exiting ...\n");
        return 1;
    }

    for (int i=0;i<3;i++) {
        c = is.peek ();
        while (c == '$')
            is.getline (line,256);
        for (int kk=0;kk<ntin;kk++) {
            is >> tph[i][kk];
        }
        // Change sign for phase, to be negative for lags
        for (int kk=0;kk<ntin;kk++) {
            tph[i][kk] = -tph[i][kk];
        }
    }
    if ( !is || is.fail() ) {
        printf ("\nError reading data lines. Exiting ...\n");
        return 1;
    }

    return 0;
}

int hardisp (const int* idate,const int& it_size,const int& irnt,
        const double& samp,const char* filename)
{
    /*+---------------------------------------------------------------------
     *
     *  Parameters below set the buffer size for computing the tides
     *  recursively (nl), the number of harmonics used in the prediction
     *  (nt; this must also be set in the subroutine admint) and the number
     *  of harmonics read in (ntin)
     *
     *----------------------------------------------------------------------*/
    constexpr int nl   ( 600 );
    constexpr int nt   ( 342 );

    double dr ( 0.01745329252e0 );
    int irli  ( 1 );
    int it[5];

    constexpr double PI = 3.1415926535897932384626433e0;

    //  Cartwright-Tayler numbers of tides used in Scherneck lists:
    //+ M2, S2, N2, K2, K1, O1, P1, Q1, Mf, Mm, Ssa
    static int idt[][6] = {
        {2, 0, 0, 0, 0, 0},
        {2, 2,-2, 0, 0, 0},
        {2,-1, 0, 1, 0, 0},
        {2, 2, 0, 0, 0, 0},
        {1, 1, 0, 0, 0, 0},
        {1,-1, 0, 0, 0, 0},
        {1, 1,-2, 0, 0, 0},
        {1,-2, 0, 1, 0, 0},
        {0, 2, 0, 0, 0, 0},
        {0, 1, 0,-1, 0, 0},
        {0, 0, 2, 0, 0, 0}
    };

    /*+----------------------------------------------------------------------
     *
     *  Read and assign the date
     *
     *-----------------------------------------------------------------------*/
    if (it_size == 5) {
        std::copy (idate,idate+5,it);
    } else if (it_size == 6) {
        std::copy (idate,idate+5,it);
        int month ( idate[1] );
        int dom   ( idate[2] );
        it[1] = dom + iers2010::hisp::mday (it[0],month);
        it[0] = idate[0];
        it[2] = idate[3];
        it[3] = idate[4];
        it[4] = idate[5];
    } else {
        return 1;
    }

    /*+---------------------------------------------------------------------
     *  Read in amplitudes and phases, in standard "Scherneck" form, from
     *  input file
     *----------------------------------------------------------------------*/
    std::ifstream fin;
    fin.open (filename,std::ifstream::in);
    if (!fin.is_open ())
        return 1;
    if ( collect_argv (fin) ) {
        fin.close ();
        return 1;
    }
    fin.close ();
     // WARNING Neet to read better; Skip lines starting with '$',  skip or read
     //  station name
    
    /*+---------------------------------------------------------------------
     *
     *  Find amplitudes and phases for all constituents, for each of the
     *  three displacements. Note that the same frequencies are returned 
     *  each time.
     *
     *  BLQ format order is vertical, horizontal EW, horizontal NS
     *
     *----------------------------------------------------------------------*/
    double az[nt],pz[nt],f[nt],aw[nt],pw[nt],as[nt],ps[nt];
    int ntout;
    double amp[ntin], phase[ntin];
    
    for (int i=0;i<ntin;i++) {
        amp[i]   = tamp[0][i];
        phase[i] = tph[0][i];
    }
    iers2010::hisp::admint (amp,idt,phase,az,f,pz,ntin,ntout,it);

    for (int i=0;i<ntin;i++) {
        amp[i] = tamp[1][i];
        phase[i] = tph[1][i];
    }
    iers2010::hisp::admint (amp,idt,phase,aw,f,pw,ntin,ntout,it);

    for (int i=0;i<ntin;i++) {
        amp[i] = tamp[2][i];
        phase[i] = tph[2][i];
    }
    iers2010::hisp::admint (amp,idt,phase,as,f,ps,ntin,ntout,it);

    // set up for recursion, by normalizing frequencies, and converting
    // phases to radians
    double wf[nt];
    for (int i=0;i<ntout;i++) {
        pz[i] = dr * pz[i];
        ps[i] = dr * ps[i];
        pw[i] = dr * pw[i];
        f[i]  = samp * PI * f[i]/43200.e0;
        wf[i] = f[i];
    }
    
    /*+---------------------------------------------------------------------
     *
     *  Loop over times, nl output points at a time. At the start of each
     *  such block, convert from amp and phase to sin and cos (hc array) at
     *  the start of the block. The computation of values within each
     *  block is done recursively, since the times are equi-spaced.
     *
     *----------------------------------------------------------------------*/
    while ( true ) {
        int irhi ( std::min (irli+nl-1,irnt) );
        int np   ( irhi - irli + 1 );

        // Set up harmonic coefficients, compute tide, and write out
        double hcz[2*nt+1],hcs[2*nt+1],hcw[2*nt+1];
        for (int i=0;i<nt;i++) {
            hcz[2*i] = az[i] * cos (pz[i]);
            hcz[2*i+1]  = -az[i] * sin (pz[i]);
            hcs[2*i] = as[i] * cos (ps[i]);
            hcs[2*i+1]  = -as[i] * sin (ps[i]);
            hcw[2*i] = aw[i] * cos (pw[i]);
            hcw[2*i+1]  = -aw[i] * sin (pw[i]);
        }

        double dz[nl],ds[nl],dw[nl];
        double scr[3*nt];
        iers2010::hisp::recurs (dz,np,hcz,ntout,wf,scr);
        iers2010::hisp::recurs (ds,np,hcs,ntout,wf,scr);
        iers2010::hisp::recurs (dw,np,hcw,ntout,wf,scr);

        for (int i=0;i<np;i++)
            printf ("\n[%02i] %+14.6f %+14.6f %+14.6f",i+1,dz[i],ds[i],dw[i]);
        printf ("\n");

        if (irhi==irnt)
            break;

        irli = irhi + 1;

        // Reset phases to the start of the new section
        for (int i=0;i<nt;i++) {
            pz[i] = fmod (pz[i] + np * f[i],2.e0*PI);
            ps[i] = fmod (ps[i] + np * f[i],2.e0*PI);
            pw[i] = fmod (pw[i] + np * f[i],2.e0*PI);
        }
    }

    return 0;
}
