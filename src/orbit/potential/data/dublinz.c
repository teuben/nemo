/*
 * dublinz.c: potential & forces as defined by a rotation curve
 *            that is linear to (r0,v0) and flat thereafter
 *          and quasi harmonic in Z
 *
 *      06-may-03       RPO, adapted from flatz
 *      19-sep-04       PJT: double/float
 */

/*CTEX
 *  {\bf potname=dublinz
 *       potpars={\it $\Omega,r_0,r_1,v_1,dvdr,s,h$}}
 *	 
 * Forces defined by a double linear rotation curve defined by
 * ($r_1,v_1$) and a gradient $dvdr$ between $r_0$ and $r_1$.
 * As in  {\bf flatz} (from which this one is derived), the 
 * potential is quasi harmonic in $Z$ (linear forces), 
 * with radial scalelength $h$ and scale height $s$
 */
 

#include <stdinc.h>
#include <spline.h>
#include <table.h>
#include <potential_float.h>

local double omega = 0.0;
local double r0    = 0.001; /* inner radius  */
local double r1    = 1.00;  /* outer radius  */
local double v1    = 1.0;   /* rotation speed beyond r0                      */
local double dvdr  = 0.0;   /* RC gradient beyond r0                         */
local double s     = 1.0;   /* columndensity of matter (constant with r)     */
local double h     = 0.2;   /* scaleheight of the vertical mass distributio  */
local double h2;

void inipotential (int *npar, double *par, string name)
{
    int n;

    n = *npar;
    if (n>0) omega = par[0];
    if (n>1) r0    = par[1];
    if (n>2) r1    = par[2];
    if (n>3) v1    = par[3];
    if (n>4) dvdr  = par[4];
    if (n>5) s     = par[5];
    if (n>6) h     = par[6];
    if (n>7) warning("DublinZ potential: only 6 parameters usable");
    
    dprintf (1,"INIPOTENTIAL DublinZ potential %s\n",name);
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  R0 = %g R1 = %g V1 = %g dVdR= %g S = %g  H = %g\n", r0, r1, v1, dvdr, s, h);

    par[0] = omega;     /* return pattern speed again */
    h2     = h*h;
}
    
void potential_double (int *ndim, double *pos,double *acc,double *pot,double *time)
{
    real r, r2, v, v0, f;
    int    i;

    for (i=0, r2=0.0; i<2; i++)
        r2 += sqr(pos[i]);
    r=sqrt(r2);

    v0 = v1 - (r1-r0)*dvdr;

    if (r < r0)
        v = (r/r0)*v0;
    else
      v = v0 + (r-r0)*dvdr;

    if (r > 0)
	f = sqr(v/r);
    else
    	f = 0;

    *pot   =  0.0;             /* no potentials... for now */
    acc[0] = -f*pos[0]; 
    acc[1] = -f*pos[1]; 
    acc[2] = -s*pos[2]/sqrt(pos[2]*pos[2]+h2);    /* --> 1/sqrt(z^2+h^2) */
}
