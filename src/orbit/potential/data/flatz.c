/*
 * flatz.c: potential & forces as defined by a rotation curve
 *            that is linear to (r0,v0) and flat thereafter
 *          and quasi harmonic in Z
 *
 *	29-dec-01	derived from rotcur
 *      13-feb-03       derived from rotcur0; toy project Rob Olling
 *      15-feb-03       RPO, buglet fixed and documentation added
 *      19-sep-04       PJT, double/float
 */

/*CTEX
 *  {\bf potname=flatz
 *       potpars={\it $\Omega,r_0,v_0,s,h$}}
 *
 * forces defined by a rotation curve that is linear to 
 * $(r_0,v_0)$ and flat thereafter and quasi harmonic in $Z$,
 * with radial scalelength $h$ and scale height $s$.
 * See also {\bf dublinz} for a variation on this theme.
 *	 
 */
 

#include <stdinc.h>
#include <spline.h>
#include <table.h>
#include <potential_float.h>

local double omega = 0.0;
local double r0    = 1.0;  /* radius at which RC changes from linear to flat */
local double v0    = 1.0;  /* rotation speed beyond r0                       */
local double s     = 1.0;  /* columndensity of matter (constant with r)      */
local double h     = 0.2;  /* scaleheight of the vertical mass distribution  */
local double h2;

void inipotential (int *npar, double *par, string name)
{
    int n;

    n = *npar;
    if (n>0) omega = par[0];
    if (n>1) r0    = par[1];
    if (n>2) v0    = par[2];
    if (n>3) s     = par[3];
    if (n>4) h     = par[4];
    if (n>5) warning("FlatZ potential: only 5 parameters usable");
    
    dprintf (1,"INIPOTENTIAL FlatZ potential %s\n",name);
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  R0 = %g  V0 = %g S = %g  H = %g\n", r0, v0, s, h);

    par[0] = omega;     /* return pattern speed again */
    h2     = h*h;
}
    
void potential_double (int *ndim, double *pos,double *acc,double *pot,double *time)
{
    real r, r2, v, f;
    int    i;

    for (i=0, r2=0.0; i<2; i++)
        r2 += sqr(pos[i]);
    r=sqrt(r2);

    if (r < r0)
        v = (r/r0)*v0;
    else
        v = v0;

    if (r > 0)
	f = sqr(v/r);
    else
    	f = 0;

    *pot   =  0.0;             /* no potentials... for now */
    acc[0] = -f*pos[0]; 
    acc[1] = -f*pos[1]; 
    acc[2] = -s*pos[2]/sqrt(pos[2]*pos[2]+h2);    /* --> 1/sqrt(z^2+h^2) */
}
