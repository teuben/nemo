/*  
 * POTENTIAL C_TO_F
 *
 * stub in order for C to call Fortran routines.
 * only works on the the BSD convention of F2C interface
 *
 */


#include <stdinc.h>

void inipotential(int *npar, double *par, string name)
{
    inipotential_(npar, par, name, strlen(name));
}

void potential (int *ndim,double *pos,double *acc,double *pot,double *time)
{
    potential_(ndim,pos,acc,pot,time);
}
