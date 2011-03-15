/* 
 * CCDPRINT: print values off gridpoints
 *
 *	 15-mar-2011   V1.0 created        Peter Teuben
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
	"in=???\n  Input filename",
	"x=0\n              Pixels in X to print (0=1st pixel)",
	"y=0\n              Pixels in Y to print",
	"z=0\n              Pixels in Z to print",
	"scale=1.0\n        Scale factor for printout",
	"format=%g\n        Format specification for output",
	"newline=f\n        Force newline between each line?",
	"label=\n	    Add x, y and or z labels add appropriate labels",
	"offset=0\n         Offset (0 or 1) to index coordinates X,Y,Z",
	"pixel=f\n          XYZ in Pixel or Physical coordinates?",
	"mode=2\n           2d or 3d interpolation?",
	"VERSION=1.0\n      15-mar-2011 PJT",
	NULL,
};

string usage = "interpolate values from an image";

string cvsid="$Id$";

int ini_array(string key, real *dat, int ndat, real offset);
void myprintf(string fmt, real v);


nemo_main()
{
    int     i, j, k, j1, l, ix,iy,iz,nx, ny, nz, nxpos, nypos, nzpos, nlcount=0;
    real    *ixr, *iyr, *izr, offset, dx, dy, dz;
    int     nmax, mode;
    bool    newline, xlabel, ylabel, zlabel, Qpixel;
    string  infile;			        /* file name */
    stream  instr;				/* file stream */
    imageptr iptr=NULL;			      /* allocated dynamically */
    string   fmt, label;
    real     scale_factor, x, y, z, f, f00, f01, f10, f11;

    scale_factor = getdparam("scale");
    fmt = getparam("format");
    instr = stropen (getparam("in"), "r");
    newline = getbparam("newline");
    Qpixel = getbparam("pixel");
    label = getparam("label");
    xlabel = scanopt(label,"x");
    ylabel = scanopt(label,"y");
    zlabel = scanopt(label,"z");
    offset = getrparam("offset");
    mode = getiparam("mode");
    if (read_image (instr,&iptr) == 0)
      error("Problem reading image from in=",getparam("in"));
    if(mode != 2) error("3d not implemented");
    
    nx = Nx(iptr);	                        /* cube dimensions */
    ny = Ny(iptr);
    nz = Nz(iptr);
    ixr = (real *) allocate(nx*sizeof(real));        /* allocate position arrays */
    iyr = (real *) allocate(ny*sizeof(real));
    izr = (real *) allocate(nz*sizeof(real));
    nxpos = ini_array("x",ixr,nx,offset);
    nypos = ini_array("y",iyr,ny,offset);
    nzpos = ini_array("z",izr,nz,offset);
    if (nxpos != nypos || nypos != nzpos) error("Not same number of x,y,z");
    for (l=0; l<nxpos; l++) {
      if (Qpixel) {
	ix = ixr[l];   
	iy = iyr[l];   
	iz = izr[l];   
      } else {
	ix = (ixr[l] - Xmin(iptr))/Dx(iptr);
	iy = (iyr[l] - Ymin(iptr))/Dy(iptr);
	iz = (izr[l] - Zmin(iptr))/Dz(iptr);
      }
      if (ix<0 || ix>nx-1) error("x=%g out of bounds",ixr[l]);
      if (iy<0 || iy>ny-1) error("y=%g out of bounds",iyr[l]);
      if (iz<0 || iz>nz-1) error("z=%g out of bounds",izr[l]);
      dx = ixr[l] - ix;
      dy = iyr[l] - iy;
      dz = izr[l] - iz;

      f00 = CubeValue(iptr,ix,iy,iz);
      f01 = CubeValue(iptr,ix,iy+1,iz);
      f10 = CubeValue(iptr,ix+1,iy,iz);
      f11 = CubeValue(iptr,ix+1,iy+1,iz);

      f = f00;  /* plus a bit to 2d interpolate from */

      if (Qpixel) {
	x = Xmin(iptr) + ix * Dx(iptr);
	y = Ymin(iptr) + iy * Dy(iptr);
	z = Zmin(iptr) + iz * Dz(iptr);
      } else {
	x = ixr[l];
	y = iyr[l];
	z = izr[l];
      }
      if (xlabel) myprintf(fmt,x);
      if (ylabel) myprintf(fmt,y);
      if (zlabel) myprintf(fmt,z);
      printf (fmt,f*scale_factor);	  
      printf(" ");
      if (newline) printf("\n");
    }
    if (!newline) printf("\n");
    strclose(instr);
}

ini_array(
	  string key,             /* keyword */
	  real *dat,              /* array */
	  int ndat,               /* length of array */
	  real offset)            /* offset applied to index array */
{
    int i, n;

    if (ndat <= 0) return ndat;
    n = hasvalue(key) ? nemoinpr(getparam(key),dat,ndat) : 0;
    if (n > 0) {
        for (i=0; i<n; i++) {
            dat[i] -= offset;
            if (dat[i] < 0 || dat[i] >= ndat)
	      error("Illegal boundary in %s",key);
        }
    } else
      error("Error %d parsing %s=%s",n,key,getparam(key));
    return n;
}

void myprintf(string fmt,real v)
{
    printf(fmt,v);
    printf(" ");
}
