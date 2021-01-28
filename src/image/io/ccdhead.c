/* 
 *	CCDHEAD: print out image header
 *
 *      27-apr-2011   Created , finally
 *      27-jan-2021   fixed for axis=1
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <extstring.h>
#include <image.h>

string defv[] = {
  "in=???\n       Input image filename",
  "VERSION=1.3\n  27-jan-2021 PJT",
  NULL,
};

string usage="print out image header";

string cvsid="$Id$";

void aminmax(real xm, real xr, real dx, int n, real *xmin, real *xmax) {
  *xmin = xm - (xr+0.5)*dx;
  *xmax = xm + (n-0.5-xr)*dx;
}

void nemo_main()
{
    imageptr iptr=NULL, optr=NULL;
    string *filename;
    stream instr, outstr;
    int i, j, i0, j0, nfiles;
    int nx, ny, nx1, ny1, ix, iy, n;
    real xmin,xmax,ymin,ymax,zmin,zmax;

    filename = getparam("in");
    instr = stropen(filename,"r");
    read_image(instr,&iptr);
    strclose(instr); 

    aminmax(Xmin(iptr),Xref(iptr),Dx(iptr),Nx(iptr),  &xmin, &xmax);
    aminmax(Ymin(iptr),Yref(iptr),Dy(iptr),Ny(iptr),  &ymin, &ymax);
    aminmax(Zmin(iptr),Zref(iptr),Dz(iptr),Nz(iptr),  &zmin, &zmax);

    printf("Size:      %d %d %d\n", Nx(iptr), Ny(iptr), Nz(iptr));
    printf("Cell:      %g %g %g\n", Dx(iptr), Dy(iptr), Dz(iptr));
    if (Axis(iptr) == 0) {
      printf("AXIS=0:\n");
      printf("LL-Corner: %g %g %g\n", Xmin(iptr), Ymin(iptr), Zmin(iptr));
      printf("TR-Corner: %g %g %g\n", Xmin(iptr)+(Nx(iptr)-1)*Dx(iptr),
	     Ymin(iptr)+(Ny(iptr)-1)*Dy(iptr),
	     Zmin(iptr)+(Nz(iptr)-1)*Dz(iptr));
      printf("X-range:   %g %g\n", xmin,xmax);
      printf("Y-range:   %g %g\n", ymin,ymax);
      printf("Z-range:   %g %g\n", zmin,zmax);
    } else {
      printf("AXIS=1:\n");
      printf("LL-Corner: %g %g %g\n",
	     Xmin(iptr)-(Nx(iptr)-1)*Dx(iptr)/2,
	     Ymin(iptr)-(Ny(iptr)-1)*Dy(iptr)/2,
	     Zmin(iptr)-(Nz(iptr)-1)*Dz(iptr)/2);
	     
      printf("TR-Corner: %g %g %g\n",
	     Xmin(iptr)+(Nx(iptr)-1)*Dx(iptr)/2,
	     Ymin(iptr)+(Ny(iptr)-1)*Dy(iptr)/2,
	     Zmin(iptr)+(Nz(iptr)-1)*Dz(iptr)/2);
      printf("X-range:   %g %g\n", xmin,xmax);
      printf("Y-range:   %g %g\n", ymin,ymax);
    }
    printf("MinMax:    %g %g\n", MapMin(iptr), MapMax(iptr));
}


