/* 
 * CCD: potential (and interpolated forces) from an image : 2D version
 *      z-forces are merely returned as 0.0
 *
 *	14-oct-91  V1.0  -  created - but never finished        PJT
 *	19-oct-91  ... finished. Anything beyond an edge has no forces
 *		   working on them, and will hence move away (escape)   PJT
 *	   oct-93  get_pattern
 *	 4-mar-96  fixed VERSION to CCD_VERSION
 *      19-jul-02  cleaned up the code a bit, provided scale factor     PJT
 *
 *  Note: astronomical images often have Dx<0 (Right Ascension increases
 *        to the left)....
 *
 */

/*CTEX
 *  {\bf potname=ccd
 *	 potpars={\it $\Omega,I_scale,X_o,Y_o$}
 *	 potfile={\it image(5NEMO)}}
 */

#include <stdinc.h>
#include <filestruct.h>
#include <image.h>

#define CCD_VERSION "ccd V2.0 17-jul-02"

local double   omega = 0.0;
local double   iscale = 1.0;
local stream   potstr = NULL;
local imageptr iptr = NULL;
local double   idx, idy, xmin, ymin, xmax, ymax, dx, dy;
local int      nx, ny;
local double   xoff=0.0;       /* central pixel shift w.r.t. Xmin.Ymin */
local double   yoff=0.0;       /* in pixel coordinates */


void inipotential (int *npar, double par[], char *name)
{
    int n;

    n = *npar;
    if (n>0)  omega = par[0];
    if (n>1)  iscale = par[1];
    if (n>2)  xoff = par[2];
    if (n>3)  yoff = par[3];
    if (n>4)  warning("inipotential(ccd): npar=%d only 4 parameter accepted",n);

    dprintf(1,"INIPOTENTIAL: %s: %s\n",CCD_VERSION,name);
    dprintf(1,"  Parameters:  Omega=%g iscale=%g xoff,yoff=%g,%g\n",
	    omega,iscale,xoff,yoff);

    potstr = stropen (name, "r");          /* open the image */
    read_image (potstr,&iptr);              /* read the image */
    if (iscale != 1.0) {
      int i,j;
      for (j=0; j<Ny(iptr); j++)
      for (i=0; i<Nx(iptr); i++)
	MapValue(iptr,i,j) =  MapValue(iptr,i,j) * iscale;
    }

    dx = Dx(iptr);
    dy = Dy(iptr);
    idx = 1.0/dx;
    idy = 1.0/dy;
    xmin = Xmin(iptr);
    ymin = Ymin(iptr);
    nx = Nx(iptr);
    ny = Ny(iptr);

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
    xmax = xmin + nx * dx;
    ymax = ymin + ny * dy;
    dprintf(1,"Offset and scale factors: xmin,ymin,1/dx,1/dy=%f %f %f %f\n",
            xmin,ymin,idx,idy);
    dprintf(1,"Formal X-range: %g %g\n",xmin - dx/2,xmax + dx/2);
    dprintf(1,"Formal Y-range: %g %g\n",ymin - dy/2,ymax + dy/2);

    par[0] = omega;
}

void potential(int *ndim,double *pos,double *acc,double *pot,double *time)
{
    double x, y, dx, x3, x4, x8, dy, y3, y4, y8;
    int ix, iy;

    x = (pos[0]-xmin)*idx + 0.5 + xoff;    /* fractional cell index   0..nx */
    y = (pos[1]-ymin)*idy + 0.5 + yoff;
    ix = (int) floor(x);            /* cell index:   0 .. nx-1 */
    iy = (int) floor(y);
    dx = x-ix-0.5;    x3=3*dx;    x4=4*dx;    x8=2*x4;
    dy = y-iy-0.5;    y3=3*dy;    y4=4*dy;    y8=2*y4;

    dprintf(5,"x,y,ix,iy=%f %f %d %d\n",x,y,ix,iy);
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
    acc[2] = 0.0;
}



