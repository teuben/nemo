/*
 *  MAIN.C : NEMO driver for Hernquist' fortran treecode program
 *
 */

#include <stdinc.h>

string defv[] = {
    "in=???\nInput file",
    "out=???\nOutput file",
    "nsteps=100\nNumber of steps to integrate",
    "nout=10\nNumber of steps between major output",
    "dtime=0.05\nIntegration time step",
    "tol=1.0\nTolerance (opening angle)",
    "eps=0.05\nSoftening length",
    "usequad=f\nUse quadrupole corrections?",
    "headline=\nAny further comments for this run",
    "VERSION=1.0\n20-jul-89",
    NULL,
};

nemo_main()
{
    pars_2_lars();
    call_2_lars();                 /* call worker routine until done */
}

