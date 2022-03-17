/*
 * POTCCD: Turn potential into an (potential|density) image;
 *	   useful for potname=ccd
 *
 *	10-jun-92	Created		        Peter Teuben
 *	30-mar-94	poisson density		pjt
 *	 2-jul-95	documented bug fix on rotating potentials 	pjt
 *      19-apr-96       improved error reporting        pjt
 *      12-jun-98       made ndim=3 the default, like potlist           pjt
 *	13-sep-01	better prototype for proc			pjt
 *       4-dec-01       also compute min/max
 *      12-sep-02       optionally output a force                       pjt
 *	22-oct-02       also allow ar,at                                pjt
 *      16-mar-2021     axis=1 proper WCS                               PJT
 *       9-jan-2022     fix Xref for axis=1                             pjt
 */

#include <stdinc.h>
#include <getparam.h>
#include <image.h>
#include <potential.h>

string defv[] = {
    "out=???\n      Output file (image)",
    "potname=???\n  Name of potential",
    "potpars=\n     Parameters for potential",
    "potfile=\n     Any optional data file associated with potential",
    "x=0\n          X-coordinate(s) to test potential at",
    "y=0\n          Y-coordinate(s) to test potential at",
    "z=0\n          Z-coordinate(s) to test potential at",
    "t=0.0\n        Time to test potential at",
    "mode=pot\n     Output pot,ax,ay,az,ar,at,den",
    "dr=\n          Differential step for (Poisson) density map",
    "omega=\n       Use this instead of any returned pattern speed",
    "ndim=3\n       Poisson map using 2D or 3D derivatives",
    "nder=1\n       1: use force der   2: use twice potential der for density",
    "VERSION=2.2a\n 9-jan-2022 PJT",
    NULL,
};

string usage = "Create potential or density image from a NEMO potential";

/*                     this code uses double[3*MAXPT] */
#ifndef MAXPT
# define MAXPT 10000
#endif

local potproc_double mypot;    /* pointer to potential calculator function */

void nemo_main(void)
{
    int    nx,ny,nz, ix,iy,iz;
    double pos[3],acc[3],pot,den,dr,time;
    double xarr[MAXPT],yarr[MAXPT],zarr[MAXPT];
    double omega, dmin=0, dmax=0;
    string mode = getparam("mode");
    int idim, ndim, maxdim = 3, first=1, idx = 0;
    int axis = 0; // new WCS system
    int nder = getiparam("nder");
    imageptr iptr;
    stream ostr;

    ostr = stropen(getparam("out"),"w");
    nx = nemoinpd(getparam("x"), xarr, MAXPT);  /* get sample arrays */
    ny = nemoinpd(getparam("y"), yarr, MAXPT);
    nz = nemoinpd(getparam("z"), zarr, MAXPT);
    dprintf(0,"Creating image %d * %d * %d\n",nx,ny,nz);
    if (nx > 1 || ny > 1 || nz > 1) {  /* check if > 1 */
        if (nx > 1) {
           if (ny>1 && ny!=nx) error("ny <> nx\n");
           if (nz>1 && nz!=nx) error("nz <> nx\n");
        } else if (ny > 1) {
           if (nz>1 && nz!=ny) error("nz <> ny\n");
        } 
    } 
    if (nx < 0 || ny < 0 || nz < 0)
        error("problem with your grid, check your x=,y=,z=");
    ndim = getiparam("ndim");
    if (hasvalue("dr")) {
        dr = getdparam("dr");
        dprintf(0,"%dD-Poisson density map with delta-r=%g and nder=%d\n",ndim,dr,nder);
    } else
        dr = -1.0;
    time = getdparam("t");
    if (streq(mode,"pot"))
      idx = 0;
    else if (streq(mode,"ax"))
      idx = 1;
    else if (streq(mode,"ay"))
      idx = 2;
    else if (streq(mode,"az"))
      idx = 3;
    else if (streq(mode,"ar"))
      idx = 4;
    else if (streq(mode,"at"))
      idx = 5;
    else if (streq(mode,"den")) {
      idx = 0;
      if (dr < 0) error("Need to supply a small positive value for dr=");
    } else
      error("bad mode=%s; allowed are: ax,ay,az,ar,at,den",mode);
    mypot = get_potential(getparam("potname"), 
			  getparam("potpars"), 
			  getparam("potfile"));
			  
    if (mypot==NULL) error("Potential could not be loaded");
    if (hasvalue("omega"))
      omega = getdparam("omega");
    else
      omega = get_pattern();
    dprintf(0,"using Omega = %g\n",omega);

    create_cube(&iptr,nx,ny,nz);
    Dx(iptr) = (nx > 1 ? xarr[1]-xarr[0] : 0.0);
    Dy(iptr) = (ny > 1 ? yarr[1]-yarr[0] : 0.0);
    Dz(iptr) = (nz > 1 ? zarr[1]-zarr[0] : 0.0);

    axis = 1;
    Axis(iptr) = axis;
    if (axis == 0) {        // deprecated
      Xmin(iptr) = xarr[0];
      Ymin(iptr) = yarr[0];
      Zmin(iptr) = zarr[0];
      Xref(iptr) = 0.0;
      Yref(iptr) = 0.0;
      Zref(iptr) = 0.0;
    } else {
      dprintf(1,"X:  %g .. %g  -> %g\n", xarr[0], xarr[nx-1], -xarr[0]/Dx(iptr));
      dprintf(1,"Y:  %g .. %g  -> %g\n", yarr[0], yarr[nx-1], -yarr[0]/Dy(iptr));
      Xmin(iptr) = 0.0;  // force the reference pixel at (0,0)
      Ymin(iptr) = 0.0;
      Zmin(iptr) = 0.0;
      Xref(iptr) = -xarr[0]/Dx(iptr);
      Yref(iptr) = -yarr[0]/Dy(iptr);
      Zref(iptr) = 0.0;
    }

    for (iz=0; iz<nz; iz++) {
      pos[2] = zarr[iz];
      for (iy=0; iy<ny; iy++) {
	pos[1] = yarr[iy];
	for (ix=0; ix<nx; ix++) {
	  pos[0] = xarr[ix];
	  (*mypot)(&maxdim,pos,acc,&pot,&time);
	  if (dr > 0.0) {                                 /* Poisson */
	    if (nder==1) {
	      den = 0.0;
	      for (idim=0; idim<ndim; idim++) {
		pos[idim] += dr;
		(*mypot)(&maxdim,pos,acc,&pot,&time);
		den -= acc[idim];
		pos[idim] -= 2*dr;
		(*mypot)(&maxdim,pos,acc,&pot,&time);
		den += acc[idim];
		pos[idim] += dr;
	      }
	    } else {
	      den = -2*ndim*pot;
	      for (idim=0; idim<ndim; idim++) {
		pos[idim] += dr;
		(*mypot)(&maxdim,pos,acc,&pot,&time);
		den += pot;
		pos[idim] -= 2*dr;
		(*mypot)(&maxdim,pos,acc,&pot,&time);
		den += pot;
		pos[idim] += dr;
	      }
	      den = den/dr;
	    }
	    pot = den/dr;
	  } else if (idx==0) {                                        /* Potential */
	    if (omega != 0.0) {
	      pot -= 0.5*sqr(omega)*(sqr(pos[0])+sqr(pos[1]));
	    }
	  } else if (idx < 4) {
	    pot = acc[idx-1];
	  } else {
	    /* 2D only */
	    real vv,vr,rr;
	    vv = acc[0]*acc[0] + acc[1]*acc[1];
	    vr = acc[0]*pos[0] + acc[1]*pos[1];
	    rr = pos[0]*pos[0] + pos[1]*pos[1];
	    if (idx == 4)
	      pot = vr/sqrt(rr);
	    else
	      pot = sqrt(vv-vr*vr/rr);
	  }
	  CubeValue(iptr,ix,iy,iz) = pot;
	  if (first) {
	    dmin = dmax = pot;
	    first = 0;
	  } else {
	    dmin = MIN(dmin, pot);
	    dmax = MAX(dmax, pot);
	  }
	}
      }
    }
    MapMin(iptr) = dmin;
    MapMax(iptr) = dmax;
    write_image(ostr, iptr);
}
