/* 
 *	CCDHEAD: print out image header
 *
 *      27-apr-2011   Created , finally
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <extstring.h>
#include <image.h>

string defv[] = {
  "in=???\n       Input filename",
  "VERSION=1.1\n  20-feb-2018 PJT",
  NULL,
};

string usage="print out image header";

string cvsid="$Id$";

extern string *burststring(string,string);

void aminmax(real x, real dx, int n, real *xmin, real *xmax) {
  *xmin = x - 0.5*dx;
  *xmax = x + (n-0.5)*dx;
}

nemo_main()
{
    imageptr iptr=NULL, optr=NULL;
    string *filename;
    stream instr, outstr;
    int i, j, i0, j0, nfiles;
    int nx, ny, nx1, ny1, ix, iy, n;
    real xmin,xmax,ymin,ymax,zmin,zmax;

    filename = getparam("in");
    instr = stropen (filename,"r");
    read_image (instr,&iptr);
    strclose(instr); 

    aminmax(Xmin(iptr),Dx(iptr),Nx(iptr),  &xmin, &xmax);
    aminmax(Ymin(iptr),Dy(iptr),Ny(iptr),  &ymin, &ymax);
    aminmax(Zmin(iptr),Dz(iptr),Nz(iptr),  &zmin, &zmax);

    printf("OLD STYLE:\n");
    printf("Size:      %d %d %d\n", Nx(iptr), Ny(iptr), Nz(iptr));
    printf("Cell:      %g %g %g\n", Dx(iptr), Dy(iptr), Dz(iptr));
    printf("LL-Corner: %g %g %g\n", Xmin(iptr), Ymin(iptr), Zmin(iptr));
    printf("TR-Corner: %g %g %g\n", Xmin(iptr), Ymin(iptr), Zmin(iptr));    //FIX THIS
    printf("X-range:   %g %g\n", xmin,xmax);
    printf("Y-range:   %g %g\n", ymin,ymax);
    printf("Z-range:   %g %g\n", zmin,zmax);
    printf("MinMax:    %g %g\n", MapMin(iptr), MapMax(iptr));

}


