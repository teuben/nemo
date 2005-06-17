/* 
 * CCD: potential (and interpolated forces) from an image : 2D version
 *      z-forces are merely returned as 0.0
 *
 *	14-oct-91  V1.0  -  created - but never finished        PJT
 *	19-oct-91  ... finished. Anything beyond an edge has no forces
 *		   working on them, and will hence move away (escape)   PJT
 *	   oct-93  get_pattern
 *	 4-mar-96  fixed VERSION to CCD_VERSION
 *      19-jul-02  V2.0 cleaned up the code a bit, provided potential scale factor     PJT
 *      13-aug-02  V3.0 reworked WCS issues due to ported FITS files from e.g. miriad's potfft
 *       6-sep-02  V3.1 some flexibility in the interpolation method                     PJT
 *
 *  Note: astronomical images often have Dx<0 (Right Ascension increases
 *        to the left)....
 * 
 *  ToDo:  potential is just given from the nearest(?) pixel, no interpolation done yet
 *
 */

/*CTEX
 *  {\bf potname=ccd
 *	 potpars={\it $\Omega,I_s,X_c,Y_c,D_x,D_y$}
 *	 potfile={\it image(5NEMO)}}
 *
 *  \smallskip
 *  This potential is defined using potential values stored on a simple cartesian grid.
 *  Using bilinear interpolation the values and derivatives are
 *  computed at any point inside the grid. Outside the grid (as defined by the
 *  WCS in the header) the potential is not defined and assumed 0.
 *  The lower left pixel of an image in NEMO is defined as (0,0), with WCS values
 *  $X_{min},Y_{min}$ derived from the header. If the ($X_c,Y_c$) parameters are used, 
 *  these are the  0-based pixel coordinates of the center pixel. If ($D_x,D_y$) are used, 
 *  these are the pixel separations.
 *  To aid astronomical images where $D_x < 0$, these are interpreted as positive.
 *  Also note that potentials are generally negative, so it is not uncommon to need
 *  $I_s = -1$. Programs such as {\it potccd} can create such a {\bf ccd} grid 
 *  potential from a regular potential.
 *
 *  \smallskip
 *  Note: Since these forces are defined only in the Z=0 plane, the Z-forces are always
 *  returned as 0.
 */

#include <stdinc.h>
#include <potential_float.h>
#include <filestruct.h>
#include <image.h>

#define CCD_VERSION "ccd V3.1a 17-jun-05"

local double   omega = 0.0;
local double   iscale = 1.0;
local stream   potstr = NULL;
local imageptr iptr = NULL;
local double   idx, idy, xmin, ymin, xmax, ymax;
local bool     Qcen,Qdel;
local int      nx, ny;
local double   xcen=0.0;       /* central pixel; where (0,0) is the lower left pixel */
local double   ycen=0.0;
local double   dx=1;           /* pixel steps */
local double   dy=1;
local int      method = 1;


void inipotential (int *npar, double par[], char *name)
{
    int n;

    n = *npar;
    if (n>0)  omega = par[0];                /* standard pattern speed */
    if (n>1)  iscale = par[1];               /* scaling factor applied to potential */
    if (n>2)  {                              /* alternate definition of center pixel */
      xcen = par[2];
      if (n>3)  ycen = par[3];
      else ycen = xcen;
    }
    if (n>4)  {                              /* alternate definition of pixel size */
      dx = par[4];
      if (n>5)  dy = par[5];
      else dy = dx;
    }
    if (n>6)  warning("inipotential(ccd): npar=%d only 6 parameter accepted",n);

    /* set some easy to use booleans */
    Qcen = n>2;     /* if a new center was defined via parameters */
    Qdel = n>4;     /* if a new pixel size was defined via parameters */

    dprintf(1,"INIPOTENTIAL: %s: %s\n",CCD_VERSION,name);
    dprintf(1,"  Parameters:  Omega=%g iscale=%g xcen,ycen=%g,%g dx,dy=%g,%g\n",
	    omega,iscale,xcen,ycen,dx,dy);

    potstr = stropen (name, "r");          /* open the image */
    read_image (potstr,&iptr);              /* read the image */
    if (iscale != 1.0) {
      int i,j;
      for (j=0; j<Ny(iptr); j++)
      for (i=0; i<Nx(iptr); i++)
	MapValue(iptr,i,j) =  MapValue(iptr,i,j) * iscale;
    }

    nx = Nx(iptr);
    ny = Ny(iptr);
    if (!Qdel) {
      dx = Dx(iptr);
      dy = Dy(iptr);
    }
    idx = 1.0/dx;
    idy = 1.0/dy;
    if (!Qcen) {
      xmin = Xmin(iptr);
      ymin = Ymin(iptr);
    } else {
      xmin = -xcen*dx;
      ymin = -ycen*dy;
    }

    if (idx != idy) {
        if (idx == -idy && xmin == -ymin) {     /* try and patch it */
            idx = -idx;
            xmin = -xmin;
            warning("Astronomical coordinate convention assumed");
        } else
            warning("Possible bug when using dx != dy");
    }
    if (idx<0) warning("1/Dx=%f",idx);
    if (idy<0) warning("1/Dy=%f",idy);
    xmax = xmin + nx * dx;              /* pixel centers */
    ymax = ymin + ny * dy;
    dprintf(1,"Offset and scale factors: xmin,ymin,1/dx,1/dy=%f %f %f %f\n",
            xmin,ymin,idx,idy);
    dprintf(1,"Formal full pixel X-range: %g %g\n",xmin-0.5*dx,xmax+0.5*dx);
    dprintf(1,"Formal full pixel Y-range: %g %g\n",ymin-0.5*dy,ymax+0.5*dy);
    dprintf(1,"Ranges: %g %g %g %g\n",xmin,xmax,ymin,ymax);

    par[0] = omega;
}

void potential_double(int *ndim,double *pos,double *acc,double *pot,double *time)
{
    double x, y, dx, x3, x4, x8, dy, y3, y4, y8;
    int ix, iy;

    x = (pos[0]-xmin)*idx + 0.5;    /* fractional cell index   0..nx */
    y = (pos[1]-ymin)*idy + 0.5;
    ix = (int) floor(x);            /* cell index:   0 .. nx-1 */
    iy = (int) floor(y);
    dx = x-ix-0.5;    x3=3*dx;    x4=4*dx;    x8=2*x4;
    dy = y-iy-0.5;    y3=3*dy;    y4=4*dy;    y8=2*y4;


#if 0
    if (ix<1) ix=1;         /* kludge for border cell particles ??*/
    if (ix>nx-2) ix=nx-2
    if (iy<1) iy=1;
    if (iy>ny-2) iy=ny-2
#else
    if (ix<1 || ix>nx-2 || iy<1 || iy>ny-2) {	/* this will make anything */
        *pot = 0.0;				/* beyond the edge escape */
        acc[0] = acc[1] = acc[2] = 0.0;
        return;
    }
#endif

    *pot = MapValue(iptr,ix,iy);               /* no interpolation for now */

	/* args.... i don't remember where i got this wonderful formulae from */
	/* or to quote from the original Unix kernel : 			      */
	/* "you are not expected to understand this"   			      */

    
    switch (method) {
    case 1:
      acc[0] = -( (-2+x4+y3)*MapValue(iptr,ix-1,iy-1)
               +(  -x8   )*MapValue(iptr,ix  ,iy-1)
               +( 2+x4-y3)*MapValue(iptr,ix+1,iy-1)
               +(-2+x4   )*MapValue(iptr,ix-1,iy)
               +(  -x8   )*MapValue(iptr,ix  ,iy)
               +( 2+x4   )*MapValue(iptr,ix+1,iy)
               +(-2+x4-y3)*MapValue(iptr,ix-1,iy+1)
               +(  -x8   )*MapValue(iptr,ix  ,iy+1)
	       +( 2+x4+y3)*MapValue(iptr,ix+1,iy+1) )/12.0*idx;

      acc[1] = -( (-2+x3+y4)*MapValue(iptr,ix-1,iy-1)
               +(-2   +y4)*MapValue(iptr,ix  ,iy-1)
               +(-2-x3+y4)*MapValue(iptr,ix+1,iy-1)
               +(     -y8)*MapValue(iptr,ix-1,iy)
               +(     -y8)*MapValue(iptr,ix  ,iy)
               +(     -y8)*MapValue(iptr,ix+1,iy)
               +( 2-x3+y4)*MapValue(iptr,ix-1,iy+1)
               +( 2   +y4)*MapValue(iptr,ix  ,iy+1)
	       +( 2+x3+y4)*MapValue(iptr,ix+1,iy+1) )/12.0*idy;
      break;
      acc[0] = 0.0;
      acc[1] = 0.0;
    case 2:
    default:
      error("%d: Unknown method of interpolation",method);
    }
    acc[2] = 0.0;
    dprintf(5,"x,y,ix,iy,ax,ay=%f %f %d %d %g %g\n",x,y,ix,iy,acc[0],acc[1]);
}



