/* 
 * VXY: potential (and interpolated forces) from an image : 2D version
 *      z-forces are merely returned as 0.0
 *
 *	10-apr-96  cloned off ccd.c
 *
 *  Note: astronomical images often have Dx<0 (Right Ascension increases
 *        to the left);
 *
 */

/*CTEX
 *  {\bf potname=flowxy
 *	 potpars={\it $\Omega,X_o,Y_o$}
 *	 potfile={\it image(5NEMO)}}
 */

#include <stdinc.h>
#include <filestruct.h>
#include <image.h>

#define CCD_VERSION "flowxy V1.0 10-apr-96"

local double   omega = 0.0;
local double   gravc = 1.0;
local stream   potstr = NULL;
local imageptr iptrx = NULL, iptry = NULL;
local double   idx, idy, xmin, ymin;
local int      nx, ny;
local double   xoff=0.0;       /* central pixel shift w.r.t. Xmin.Ymin */
local double   yoff=0.0;       /* in pixel coordinates */


void inipotential (npar, par, name)
int *npar;
double par[];
char *name;
{
    int n;

    n = *npar;
    if (n>0)  omega = par[0];
    if (n>1)  xoff = par[1];
    if (n>2)  yoff = par[2];
    if (n>3)  warning("inipotential(flowxy): npar=%d only 3 parameter accepted",n);

    dprintf(1,"INIPOTENTIAL: %s: %s\n",CCD_VERSION,name);
    dprintf(1,"  Parameters:  Omega=%g xoff,yoff=%g,%g\n",omega,xoff,yoff);

    potstr = stropen (name, "r");            /* open the image */
    read_image (potstr,&iptrx);              /* read the 1st image */
    read_image (potstr,&iptry);              /* read the 2nd image */

    idx = 1.0/Dx(iptrx);          /* set some constants to remember */
    idy = 1.0/Dy(iptrx);
    xmin = Xmin(iptrx);
    ymin = Ymin(iptrx);
    nx = Nx(iptrx);
    ny = Ny(iptrx);

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
    dprintf(1,"Offset and scale factors: xmin,ymin,1/dx,1/dy=%f %f %f %f\n",
            xmin,ymin,idx,idy);

    par[0] = omega;
}

#define VX(ix,iy)	MapValue(iptrx,ix,iy)
#define VY(ix,iy)	MapValue(iptry,ix,iy)

void potential(ndim,pos,acc,pot,time)
int    *ndim;
double pos[], acc[], *pot, *time;
{
    double x, y, dx, x3, x4, x8, dy, y3, y4, y8;
    int ix, iy;

    x = (pos[0]-xmin)*idx + 0.5 + xoff;    /* fractional cell index   0..nx */
    y = (pos[1]-ymin)*idy + 0.5 + yoff;
    ix = (int) floor(x);            /* cell index:   0 .. nx-1 */
    iy = (int) floor(y);
    dx = x-ix;			    /* index in a cell: 0.0 .. 1.0 */
    dy = y-iy;

    dprintf(5,"pos=%g %g x,y,ix,iy=%f %f %d %d %g %g\n",pos[0],pos[1],x,y,ix,iy,dx,dy);
#if 0
    if (ix<1) ix=1;         /* kludge for border cell particles ??*/
    if (ix>nx-2) ix=nx-2
    if (iy<1) iy=1;
    if (iy>ny-2) iy=ny-2
#else
    if (ix<1 || ix>nx-2 || iy<1 || iy>ny-2) {	/* this will make anything */
        *pot = -1.0;				/* beyond the edge sit still */
        acc[0] = acc[1] = acc[2] = 0.0;
        return;
    }
#endif

    *pot = -1.0;

    acc[0] = (1-dx)*(1-dy)*VX(ix,iy) + dx*(1-dy)*VX(ix+1,iy) +
              (1-dx)*dy*VX(ix,iy+1) + dx*dy*VX(ix+1,iy+1);
    		
    acc[1] = (1-dx)*(1-dy)*VY(ix,iy) + dx*(1-dy)*VY(ix+1,iy) +
              (1-dx)*dy*VY(ix,iy+1) + dx*dy*VY(ix+1,iy+1);
    		
    acc[2] = 0.0;
}



