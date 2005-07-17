/* 
 *	CCDTWICE: clone each X and Y pixels, making an image 4 times as large
 *
 *      28-may-2004   Created for Spitzer data
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>

#include <image.h>

string defv[] = {
    "in=???\n       Input filename (2d image)",
    "out=???\n      Output filename (2d image)",
    "flux=t\n       Should flux be conserved while duplicating?",
    "n=2\n          Number of times to replicated pixels",
    "VERSION=1.0\n  30-jul-04 PJT",
    NULL,
};

string usage="duplicate a 2D image";

void slice(imageptr i, imageptr o,  int n, real factor);

nemo_main()
{
    imageptr iptr=NULL, optr=NULL;
    stream instr, outstr;
    int i, j, i0, j0;
    int nx, ny;
    int n = getiparam("n");
    real factor;
    bool Qflux = getbparam("flux");

    if (n<1) error("n=%d illegal",n);

    instr = stropen (getparam("in"),"r");	/* get stream */
    read_image (instr,&iptr);               /* read image */
    strclose(instr);                        /* close image file */

    nx = n*Nx(iptr);
    ny = n*Ny(iptr);
    factor = Qflux ? 1.0/(n*n) : 1.0;

    create_image(&optr,nx,ny);
    slice(iptr, optr, n, factor);

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

void slice(imageptr i, imageptr o, int n, real factor)
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
