/*
 * potential_float:    wrapper for potential_double
 *
 *  this routine is a small cheat, and by including this file,
 *  "automatically" converts an old-style (before V5.4) potential()
 *  to the new style potential_double() / potential_float() 
 * 
 *  Performance:  testing 100,000 evaluations on potname=hom
 *                in double mode 1.20", double+float mode 1.25
 *                thus a penalty of about 4%
 *
 *  19-sep-2004	  implemented, but only for C (Fortran portability?)
 */

extern void potential_double(int *,double *,double *,double *,double *);

void potential_float (int *ndim,float *pos,float *acc,float *pot,float *time)
{
    double dpos[3], dacc[3], dpot, dtime;
    int i;

    for (i=0; i<*ndim; i++) dpos[i] = pos[i];
    dtime = *time;
    potential_double(ndim,dpos,dacc,&dpot,&dtime);
    for (i=0; i<*ndim; i++) acc[i] = dacc[i];
    *pot = dpot;
}

