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
 *	
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
    "dr=\n          Differential step for (Poisson) density map",
    "ndim=3\n       Poisson map using 2D or 3D derivatives",
    "VERSION=1.2a\n 13-sep-01 PJT",
    NULL,
};

string usage = "Create potential or density image from a NEMO potential";

#ifndef MAXPT
# define MAXPT 2048
#endif

local potproc_double mypot;    /* pointer to potential calculator function */

void nemo_main(void)
{
    int    i, nx,ny,nz, ix,iy,iz, stepx, stepy, stepz, nsteps;
    double pos[3],acc[3],pot,den,dr,da[3],time;
    double xarr[MAXPT],yarr[MAXPT],zarr[MAXPT];
    double ax,ay,az,epot;
    double fourpi = FOUR_PI;
    double omega;
    char *fmt, s[20], pfmt[256];
    imageptr iptr;
    stream ostr;
    int idim, ndim, maxdim = 3;

    ostr = stropen(getparam("out"),"w");
    nx = nemoinpd(getparam("x"), xarr, MAXPT);  /* get sample arrays */
    ny = nemoinpd(getparam("y"), yarr, MAXPT);
    nz = nemoinpd(getparam("z"), zarr, MAXPT);
    dprintf(0,"Creating image %d * %d * %d\n",nx,ny,nz);
    if (nx > 1 || ny > 1 || nz > 1) {  /* check if > 1 */
        if (nx > 1) {
           nsteps = nx;
           if (ny>1 && ny!=nx) error("ny <> nx\n");
           if (nz>1 && nz!=nx) error("nz <> nx\n");
        } else if (ny > 1) {
           nsteps = ny;
           if (nz>1 && nz!=ny) error("nz <> ny\n");
        } else
           nsteps = nz;
    } else                                 /* only one position */
        nsteps = 1;
    if (nx < 0 || ny < 0 || nz < 0)
        error("problem with your grid, check your x=,y=,z=");
    ndim = getiparam("ndim");
    if (hasvalue("dr")) {
        dr = getdparam("dr");
        dprintf(0,"%dD-Poisson density map with delta-r=%g\n",ndim,dr);
    } else
        dr = -1.0;
    time = getdparam("t");

    mypot = get_potential(getparam("potname"), 
			  getparam("potpars"), 
			  getparam("potfile"));
			  
    if (mypot==NULL) error("Potential could not be loaded");
    omega = get_pattern();
#if 1
    dprintf(0,"Omega = %g\n",omega);
#else
    dprintf(0,"Omega = %g ** ignored **\n",omega);
#endif

    create_cube(&iptr,nx,ny,nz);
    Dx(iptr) = (nx > 1 ? xarr[1]-xarr[0] : 0.0);
    Dy(iptr) = (ny > 1 ? yarr[1]-yarr[0] : 0.0);
    Dz(iptr) = (nz > 1 ? zarr[1]-zarr[0] : 0.0);
    Xmin(iptr) = xarr[0];
    Ymin(iptr) = yarr[0];
    Zmin(iptr) = zarr[0];


    for (iz=0; iz<nz; iz++) {
        pos[2] = zarr[iz];
        for (iy=0; iy<ny; iy++) {
            pos[1] = yarr[iy];
            for (ix=0; ix<nx; ix++) {
                pos[0] = xarr[ix];
                (*mypot)(&maxdim,pos,acc,&pot,&time);
                if (dr > 0.0) {                                 /* Poisson */

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
                    CubeValue(iptr,ix,iy,iz) = den/dr;
                } else {                                        /* Potential */
		    if (omega != 0.0) {
            	        pot -= 0.5*sqr(omega)*
			    (sqr(pos[0])+sqr(pos[1]));
                    }
                    CubeValue(iptr,ix,iy,iz) = pot;
                }
            }
        }
    }
    write_image(ostr, iptr);
}
