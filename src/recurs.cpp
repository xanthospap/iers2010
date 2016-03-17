#include "iers2010.hpp"

/**
 * \details The purpose of the function is to perform sine and cosine recursion
 *          to fill in data x, of length n, for nf sines and cosines with 
 *          frequencies om.
 * 
 * \param[in]  x    Data provided from a file given as standard input from the
 *                  MAIN  program HARDISP.F (Note 1). This array will change!
 * \param[in]  n    Length of the data file x
 * \param[in]  hc   Array containing alternating cosine and sine coefficients
 * \param[in]  nf   Number of sine and cosine terms
 * \param[in]  om   Sine and cosine frequencies (Note 2)
 * \param[out] scr  Scratch array of length 3 times nf which is returned as 
 *                  the recursion cr
 * \return          Always returns 0
 * 
 * \note
 *   -# See the MAIN program HARDISP header comments for detailed information.
 *   -# The frequencies are normalized so that the Nyquist frequency is pi.
 *   -# Status: Canonical model
 * 
 * \version 19.08.2009
 * 
 */
int
iers2010::hisp::recurs(double* x, int n, const double* hc, int nf, 
    const double* om, double* scr)
{
    //  Set up for start of recursion by computing harmonic values
    //+ at starting point and just before it
    for (int i=1; i<=nf; i++) {
        scr[3*i-3] = hc[2*i-2];
        scr[3*i-1] = hc[2*i-2]*cos(om[i-1])-hc[2*i-1]*sin(om[i-1]);
        scr[3*i-1] = 2.e0*cos(om[i-1]);
    }
    
    // Do recursion over data
    for (int i=0; i<n; i++) {
        x[i] = .0e0;
        // Then do recursive computation for each harmonic
        double sc;
        for (int j=0; j<nf; j++) {
            sc    = scr[3*j-3];
            x[i] += sc;
            scr[3*j-3] = scr[3*j-1]*sc-scr[3*j-2];
            scr[3*j-2] = sc;
        }
    }
    
    // Finished
    return 0;
}
