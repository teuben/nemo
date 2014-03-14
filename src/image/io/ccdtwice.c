/* 
 *	CCDTWICE: clone each X and Y pixels, making an image 4 times as large
 *
 *      28-may-2004   Created for Spitzer data
 *      20-jul-2012   prepare for fancier duplication/interpolation doubling
 *
 *  for given axlength 'm' and duplication factor 'n':
 *  in duplication mode (old V1.x) nothing complicated
 *  otherwise two options:
 *     n*m         : this will copy old cell centers, but WCS edge is now new (smaller)
 *     n*m - (n-1) : this will maintain coordinate system edges, but new pixel centers
 *
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>

#include <image.h>

string defv[] = {
    "in=???\n       Input filename (2d image)",
    "out=???\n      Output filename (2d image)",
    "flux=t\n       Should flux be conserved while duplicating?",
    "n=2\n          Number of times to replicate pixels",
    "dup=t\n        duplicate or interpolate?",
    "wcs=t\n        Retain WCS borders?",
    "ndim=2\n       2D or 3D duplication",
    "VERSION=2.0\n  20-jul-2012 PJT",
    NULL,
};

string usage="duplicate a 2D image";

void slice2(imageptr i, imageptr o,  int n, real factor, bool Qwcs);

nemo_main()
{
    imageptr iptr=NULL, optr=NULL;
    stream instr, outstr;
    int i, j, i0, j0;
    int nx, ny;
    int n = getiparam("n");
    int ndim = getiparam("ndim");
    real factor;
    bool Qflux = getbparam("flux");
    bool Qdup = getbparam("dup");
    bool Qwcs = getbparam("wcs");

    if (n<1) error("n=%d illegal",n);
    if (ndim != 2) error("ndim=%d not supported yet",ndim);

    if (!Qdup && n != 2) error("Cannot interpolate");

    instr = stropen (getparam("in"),"r");	/* get stream */
    read_image (instr,&iptr);               /* read image */
    strclose(instr);                        /* close image file */

    nx = n*Nx(iptr);
    ny = n*Ny(iptr);
    if (Qdup)
      factor = Qflux ? 1.0/(n*n) : 1.0;
    else
      factor = 1.0;   /* but figure out the factor later */

    create_image(&optr,nx,ny);
    slice2(iptr, optr, n, factor, Qwcs);

    for (j=0; j<ny; j++) {
      j0 = j/n;
      for (i=0; i<nx; i++) {
	i0 = i/n;
	MapValue(optr,i,j) = factor*MapValue(iptr,i0,j0);
      }
    }

    outstr = stropen(getparam("out"),"w");
    write_image(outstr,optr);
    strclose(outstr);
}

void slice2(imageptr i, imageptr o, int n, real factor, bool Qwcs)
{
    int x, y, z, iz;
    real fac = 1.0/n;

    Namex(o) = Namey(i);
    Namey(o) = Namez(i);
    if (Qwcs) {
      Nx(o)   = n*Nx(i);
      Ny(o)   = n*Ny(i);
      Xmin(o) = Xmin(i);
      Ymin(o) = Ymin(i);
      Dx(o) = Dx(i)*fac;
      Dy(o) = Dy(i)*fac;
    } else {
      Nx(o)   = n*Nx(i) - (n-1);
      Ny(o)   = n*Ny(i) - (n-1);
      Xmin(o) = Xmin(i) + Dx(i)*fac*(n-1)*0.5;
      Ymin(o) = Ymin(i) + Dy(i)*fac*(n-1)*0.5;
      Dx(o) = Dx(i)*fac;
      Dy(o) = Dy(i)*fac;
    }
    MapMin(o) = MapMin(i) * factor;
    MapMax(o) = MapMax(i) * factor;
}

void slice3(imageptr i, imageptr o, int n, real factor)
{
    int x, y, z, iz;
    real fac = 1.0/n;

    Namex(o) = Namey(i);
    Namey(o) = Namez(i);
    Nx(o)   = n*Nx(i);
    Ny(o)   = n*Ny(i);
    Xmin(o) = Xmin(i);
    Ymin(o) = Ymin(i);
    Dx(o) = Dx(i)*fac;
    Dy(o) = Dy(i)*fac;
    MapMin(o) = MapMin(i) * factor;
    MapMax(o) = MapMax(i) * factor;
}
