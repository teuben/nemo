/*
 * POTLIST: show some properties of a potential at a certain point in space 
 *
 *	Peter Teuben		xx-xxx-87	created
 *			V1.1	20-oct-88
 *			V2.0     9-feb-90 	added time parameter to potential()
 *						standardized keyword names pot*
 *                      V2.1    26-jun-90       more stuff for standard tables
 *			V2.2 	14-jun-91	rotating patterns
 *			V2.2a    7-mar-92       happy gcc2.0
 *                          b    2-oct-93       MAXPT a bit larger
 *                      V3.0    11-oct-93       support new get_pattern
 *			 3.1    15              omega= to allow override
 *                       3.2    17-feb-94       also allow 2D poisson
 *			    b   30-mar-94	fixed NDIM + ansi
 *			    c   22-mar-95       print out pattern speed
 *			    d   20-feb-97	adapted for SINGLEPREC
 *                          e   12-jun-98       fixed bug in NDIM=3
 */

#include <stdinc.h>
#include <getparam.h>
#include <potential.h>
#include <vectmath.h>

string defv[] = {
    "potname=???\n  Name of potential",
    "potpars=\n     Parambeters for potential (1st one is pattern speed)",
    "potfile=\n     Any optional data file associated with potential",
    "x=0\n          X-coordinate(s) to test potential at",
    "y=0\n          Y-coordinate(s) to test potential at",
    "z=0\n          Z-coordinate(s) to test potential at",
    "t=0.0\n        Time to test potential at",
    "dr=\n          Differential step to use if to test Poisson Eq.",
    "omega=\n       Use this instead of any returned pattern speed",
    "format=%g\n    Format used to print numbers",
    "ndim=3\n       Poission test in 3-dim  (XYZ) or 2-dim (XY)",
    "VERSION=3.2e\n 12-jun-98 PJT",
    NULL,
};

string usage = "show a NEMO potential";

#ifndef MAXPT
#define MAXPT 10001
#endif

proc mypot;             /* pointer to potential calculator function */

void nemo_main(void)
{
    int    i, ndim, nx,ny,nz, ix,iy,iz, stepx, stepy, stepz, nsteps;
    double pos[3],acc[3],pot,dr,da[3],time;
    real   xarr[MAXPT],yarr[MAXPT],zarr[MAXPT];
    double ax,ay,az,epot;
    double fourpi = 4*PI;
    double omega;
    bool Qdens;
    char *fmt, s[20], pfmt[256];

    Qdens = hasvalue("dr");
    if (Qdens) dr = getdparam("dr");

    nx = nemoinpr(getparam("x"), xarr, MAXPT);  /* get arrays */
    ny = nemoinpr(getparam("y"), yarr, MAXPT);
    nz = nemoinpr(getparam("z"), zarr, MAXPT);
    if (nx<0 || ny<0 || nz<0)
	warning("Error decoding x,y,z, %d=too many? (%d,%d,%d)",MAXPT,nx,ny,nz);
    stepx = (nx > 1 ? 1 : 0);                   /* get steps */
    stepy = (ny > 1 ? 1 : 0);
    stepz = (nz > 1 ? 1 : 0);
    if (nx > 1 || ny > 1 || nz > 1) {  /* check if > 1 */
        if (nx > 1) {
           nsteps = nx;
           if (ny>1 && ny!=nx) error("ny <> nx (%d,%d)",ny,nx);
           if (nz>1 && nz!=nx) error("nz <> nx (%d,%d)",nz,nx);
        } else if (ny > 1) {
           nsteps = ny;
           if (nz>1 && nz!=ny) error("nz <> ny (%d,%d)",nz,ny);
        } else
           nsteps = nz;
    } else                                 /* only one position */
        nsteps = 1;

    time = getdparam("t");

    mypot = get_potential(getparam("potname"), 
			  getparam("potpars"), 
			  getparam("potfile"));
    if (mypot==NULL)
        error("Potential could not be loaded");
    if (hasvalue("omega"))
        omega = getdparam("omega");
    else
        omega=get_pattern();
    dprintf(1,"Found pattern speed = %g\n",omega);

    fmt = getparam("format");
    strcpy(s,fmt);    /* use format from command line */
    if (strchr(fmt,' ')==NULL && strchr(fmt,',')==NULL)
        strcat(s," ");      /* append separator if none specified */
    if (Qdens) {
      sprintf(pfmt,"%s%s%s%s%s%s%s%s%s%s%s%s%s\n",s,s,s,s,s,s,s,s,s,s,s,s,s);
      dprintf(0,"x y z ax ay az phi phixx phiyy phizz rho dr time\n");
    } else {
      sprintf(pfmt,"%s%s%s%s%s%s%s%s\n",s,s,s,s,s,s,s,s);
      dprintf(0,"x y z ax ay az phi time\n");
    }
    ndim = getiparam("ndim");
    if (ndim != 3 && ndim != 2) error("NDIM=%d must be 2 or 3",ndim);
                  
    for (i=0,ix=0,iy=0,iz=0;i<nsteps;i++) {
        pos[0] = xarr[ix];
        pos[1] = yarr[iy];
        pos[2] = zarr[iz];  /* formally not used for ndim=2 */

        (*mypot) (&ndim,pos,acc,&pot,&time);	/* get forces and potential */
        ax = acc[0]; ay = acc[1]; az = acc[2];	        /* save forces */
        epot = pot;					/* and potential */
	epot -= 0.5*sqr(omega)*
			(sqr(pos[0])+sqr(pos[1])+sqr(pos[2]));

        if (Qdens) {
            da[0]=acc[0]; da[1]=acc[1]; da[2]=acc[2];	/* store forces */

            pos[0] = xarr[ix]+dr;
            (*mypot) (&ndim,pos,acc,&pot,&time);
            da[0] -= acc[0];			/* force derivative along x */
  
            pos[0]=xarr[ix]; 
            pos[1]=yarr[iy]+dr;
            (*mypot) (&ndim,pos,acc,&pot,&time);
            da[1] -= acc[1];			/* force derivative along y */

            pos[1]=yarr[iy]; 
            pos[2]=zarr[iy]+dr;
            (*mypot) (&ndim,pos,acc,&pot,&time);
            if (ndim==3) da[2] -= acc[2];	/* force derivative along z */
            else da[2] = 0.0;

            da[0] /= dr; da[1] /= dr; da[2] /= dr;
        
            printf (pfmt,
                xarr[ix], yarr[iy], zarr[iy],		/* positions */
                ax,       ay,       az,        		/* forces */
                epot,					/* potential */
                da[0], da[1], da[2],                    /* higher derivs */
                (da[0]+da[1]+da[2])/fourpi,             /* poissonian density */
                dr,      				/* difference */
                time);                                  /* time */
        } else {
            printf (pfmt,
                xarr[ix], yarr[iy], zarr[iy],		/* positions */
                ax,       ay,       az,        		/* forces */
                epot,					/* potential */
                time);                                  /* time */
        } 
        ix += stepx; iy += stepy; iz += stepz;
    }
}
