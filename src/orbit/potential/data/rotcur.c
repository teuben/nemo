/*
 * rotcur.c: potential & forces as defined by a rotation curve
 *           This version can only compute one version; i.e.
 *           on re-entry of inipotential(), old versions are lost
 *
 *	   nov-90  created?	PJT
 *      19-jul-92  finished for Bikram     PJT
 *	   oct-93  get_pattern
 *	15-apr-98  oops, made forces now negative
 *      23-jun-01  gcc warnings
 *      26-feb-03  potential now contains the "rotation curve" (or whatever)
 */

/*CTEX
 *  {\bf potname=rotcur
 *       potpars={\it $\Omega$}
 *	 potfile={\it table(5NEMO)}}
 *
 * The forces returned are the axisymmetric forces as defined by
 * a rotation curve as defined by a table given from an ascii table.
 * The potential is not computed, instead the interpolated rotation
 * curve is returned in as the potential value.
 *
 * This version can only compute one version; i.e.
 * on re-entry of inipotential(), old versions are lost.
 *  
 */
 

#include <stdinc.h>
#include <spline.h>
#include <table.h>
#include <potential_float.h>

extern int nemo_file_lines(string,int);
 
local double omega = 0.0;           /* just put to zero until implemented */

local real *rad, *vel, *coef;
local int nrad, nmax;
local int entries=0;

void inipotential (int *npar, double *par, string name)
{
    int i, n, colnr[2];
    real *coldat[2];
    stream instr;
    real rscale=1.0, vscale=1.0;
    int rcol=1, vcol=2;

    n = *npar;
    if (n>0) omega = par[0];
    if (n>1) vscale = par[1];
    if (n>2) rscale = par[2];
    if (n>3) vcol = (int) par[3];
    if (n>4) rcol = (int) par[4];
    if (n>5) warning("Rotcur potential: only 5 parameters usable: omega,vscale,rscale,vcol,rcol");
    
    if (entries>0) {
        warning("Re-entering rotcur potential(5NEMO): removed previous tables");
        free(rad);
        free(vel);
        free(coef);
    }
    entries++;

    dprintf (1,"INIPOTENTIAL Rotcur potential %s\n",name);
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  Table = %s\n",name);
    dprintf (1,"  Vscale = %g  Rscale = %g\n", vscale, rscale);
    dprintf (1,"  Vcol = %d  Rcol = %d\n", vcol, rcol);

    nmax = nemo_file_lines(name,0);
    if (nmax<=0) error("file_lines returned %d lines in %s\n",
             nmax,name);
    dprintf (1,"  Nmax = %d\n",nmax);
    rad = (real *) allocate(nmax * sizeof(real));
    vel = (real *) allocate(nmax * sizeof(real));
    coldat[0] = rad;        colnr[0] = rcol;
    coldat[1] = vel;        colnr[1] = vcol;

    instr = stropen(name,"r");
    nrad = get_atable(instr,2,colnr,coldat,nmax);
    strclose(instr);
    if (nrad==0) error("No lines (%d)read from %s",nrad,name);
    if (nrad<0) {
	nrad = -nrad;
	warning("only read part of table");
    }
    for (i=0; i<nrad; i++) {
      rad[i] *= rscale;
      vel[i] *= vscale;
    }
    coef = (real *) allocate(nrad*3*sizeof(real));
    spline(coef, rad, vel, nrad);
    par[0] = omega;
    dprintf(2,"rotcur[1]: %g %g\n",rad[0],vel[0]);
    dprintf(2,"rotcur[%d]: %g %g\n",nrad,rad[nrad-1],vel[nrad-1]);
}
    
void potential_double (int *ndim, double *pos,double *acc,double *pot,double *time)
{
    real r, r2, v, f;
    int    i;

    for (i=0, r2=0.0; i<*ndim; i++)
        r2 += sqr(pos[i]);
    r=sqrt(r2);
    if (r<rad[0] || r>rad[nrad-1]) {
	*pot = 0.0;
	acc[0] = acc[1] = acc[2] = 0.0;
        return;
    } 

    v = seval( r, rad, vel, coef, nrad);     /* interpolate */
    if (r > 0)
	f = sqr(v/r);
    else
    	f = 0;
    dprintf(2,"r=%g v=%g f=%g\n",r,v,f);

    *pot = v;             /* can't compute potentials .. so just stuff the rotation curve */
    acc[0] = -f*pos[0]; 
    acc[1] = -f*pos[1]; 
    acc[2] = 0.0;
}
