/* 
 * CCDPRINT: interpolate values off gridpoints of an image
 *
 *	 15-mar-2011   V1.0 created        Peter Teuben
 *       28-apr-2011   V1.1 fixed location bug      PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
  "in=???\n           Input filename",
  "x=0\n              Pixels/Coords in X to print (0=1st pixel)",
  "y=0\n              Pixels/Coords in Y to print",
  "z=0\n              Pixels/Coords in Z to print",
  "scale=1.0\n        Scale factor for printout",
  "format=%g\n        Format specification for output",
  "newline=t\n        Force newline between each line?",
  "label=x,y\n	      Add x, y and or z labels add appropriate labels",
  "offset=0\n         Offset (0 or 1) to index coordinates X,Y,Z",
  "pixel=f\n          XYZ in Pixel? else Physical coordinates",
  "dim=2\n            2d or 3d interpolation? (use 0 for nearest)",
  "VERSION=1.1a\n      28-apr-2011 PJT",
  NULL,
};

string usage = "interpolate values from an image";

string cvsid="$Id$";

int ini_array(string key, real *dat, int ndat, real offset, bool Qpixel);
void myprintf(string fmt, real v);
int near_edge(int, int, int, int, int, int, int, int);

#ifndef MAXP
#define MAXP 100000
#endif

void nemo_main(void)
{
    int     i, j, k, j1, l, ix,iy,iz,nx, ny, nz, nxpos, nypos, nzpos;
    real    xr[MAXP], yr[MAXP], zr[MAXP], offset, dx, dy, dz;
    int     nmax, dim, nout;
    bool    newline, xlabel, ylabel, zlabel, Qpixel;
    string  infile;			        /* file name */
    stream  instr;				/* file stream */
    imageptr iptr=NULL;			      /* allocated dynamically */
    string   fmt, label;
    real     scale_factor, x, y, z, f, 
             f000, f010, f100, f110, 
             f001, f011, f101, f111;

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
    dim = getiparam("dim");
    if (read_image (instr,&iptr) == 0)
      error("Problem reading image from in=",getparam("in"));

    dprintf(1,"Map range in X: %g %g [%d]\n",
	    Xmin(iptr)-Dx(iptr)/2, Xmin(iptr)+(Nx(iptr) - 0.5)*Dx(iptr),Nx(iptr));
    dprintf(1,"Map range in Y: %g %g [%d]\n",
	    Ymin(iptr)-Dy(iptr)/2, Ymin(iptr)+(Ny(iptr) - 0.5)*Dy(iptr),Ny(iptr));
    dprintf(1,"Map range in Z: %g %g [%d]\n",
	    Zmin(iptr)-Dz(iptr)/2, Zmin(iptr)+(Nz(iptr) - 0.5)*Dz(iptr),Nz(iptr));
    
    nx = Nx(iptr);	                        /* cube dimensions */
    ny = Ny(iptr);
    nz = Nz(iptr);
#if 0
    ixr = (real *) allocate(nx*sizeof(real));        /* allocate position arrays */
    iyr = (real *) allocate(ny*sizeof(real));
    izr = (real *) allocate(nz*sizeof(real));
    nxpos = ini_array("x",ixr,nx,offset,Qpixel);
    nypos = ini_array("y",iyr,ny,offset,Qpixel);
    nzpos = ini_array("z",izr,nz,offset,Qpixel);
#else
    dprintf(1,"MAXP=%d\n",MAXP);
    nxpos = ini_array("x",xr,MAXP,offset,Qpixel);
    nypos = ini_array("y",yr,MAXP,offset,Qpixel);
    nzpos = ini_array("z",zr,MAXP,offset,Qpixel);
    if (nxpos>1 && nypos>1 && nzpos==1) {
      for (l=1; l<nxpos;l++) zr[l] = zr[0];
      nzpos = nxpos;
    }
#endif

    if (nxpos != nypos || nypos != nzpos) error("Not same number of x,y,z");
    nout = 0;
    for (l=0; l<nxpos; l++) {
      if (Qpixel) {
	ix = (int)rint(xr[l]);   
	iy = (int)rint(yr[l]);   
	iz = (int)rint(zr[l]);   
	dx = xr[l] - ix;
	dy = yr[l] - iy;
	dz = zr[l] - iz;
      } else {
	ix = (int)rint((xr[l] - Xmin(iptr))/Dx(iptr));
	iy = (int)rint((yr[l] - Ymin(iptr))/Dy(iptr));
	iz = (int)rint((zr[l] - Zmin(iptr))/Dz(iptr));
	dx = xr[l] - (Xmin(iptr)+ix*Dx(iptr));
	dy = yr[l] - (Ymin(iptr)+iy*Dy(iptr));
	dz = zr[l] - (Zmin(iptr)+iz*Dz(iptr));
      }
      if (nx==1) ix = dx = 0;
      if (ny==1) iy = dy = 0;
      if (nz==1) iz = dz = 0;

      if (near_edge(ix,iy,iz,nx,ny,nz,dim,0)) {
	nout++;
	continue;
      }

      f000 = CubeValue(iptr,ix,iy,iz);
      dprintf(2,"%d %d %d  %g %g %g %g\n",ix,iy,iz,dx,dy,dz,f000);

      if (dim == 0) {
	f = f000; 
      } else if (dim == 2) {
	if (near_edge(ix,iy,iz,nx,ny,nz,dim,1)) {
	  nout++;
	  continue;
	}
	f100 = CubeValue(iptr,ix+1,iy,iz);
	f010 = CubeValue(iptr,ix,iy+1,iz);
	f110 = CubeValue(iptr,ix+1,iy+1,iz);
	f = 
	  f000*(1-dx)*(1-dy) +
	  f100*dx*(1-dy) +
	  f010*dy*(1-dx) +
	  f110*dx*dy;
      } else if (dim==3) {
	if (near_edge(ix,iy,iz,nx,ny,nz,dim,1)) {
	  nout++;
	  continue;
	}
	f010 = CubeValue(iptr,ix,iy+1,iz);
	f100 = CubeValue(iptr,ix+1,iy,iz);
	f110 = CubeValue(iptr,ix+1,iy+1,iz);
	f001 = CubeValue(iptr,ix,iy,iz+1);
	f011 = CubeValue(iptr,ix,iy+1,iz+1);
	f101 = CubeValue(iptr,ix+1,iy,iz+1);
	f111 = CubeValue(iptr,ix+1,iy+1,iz+1);
	f = 
	  f000*(1-dx)*(1-dy)*(1-dz) +
	  f100*dx*(1-dy)*(1-dz) +
	  f010*dy*(1-dx)*(1-dz) +
	  f110*dx*dy*(1-dz);
	  f001*(1-dx)*(1-dy)*dz +
	  f101*dx*(1-dy)*dz +
	  f011*dy*(1-dx)*dz +
	  f111*dx*dy*dz;
      } else
	f = 0.0;

      if (Qpixel) {
	x = Xmin(iptr) + ix * Dx(iptr);
	y = Ymin(iptr) + iy * Dy(iptr);
	z = Zmin(iptr) + iz * Dz(iptr);
      } else {
	x = xr[l];
	y = yr[l];
	z = zr[l];
      }
      if (xlabel) myprintf(fmt,x);
      if (ylabel) myprintf(fmt,y);
      if (zlabel) myprintf(fmt,z);
      printf (fmt,f*scale_factor);	  
      printf(" ");
      if (newline) printf("\n");
    }
    if (!newline) printf("\n");
    if (nout) warning("%d point%s outside grid",nout, nout>1 ? "s" : "");
    strclose(instr);
}

int ini_array(
	  string key,             /* keyword */
	  real *dat,              /* array */
	  int ndat,               /* length of array */
	  real offset,            /* offset applied to index array */
	  bool Qpixel)            /* pixel or physical coords ? */
{
    int i, n;

    if (ndat <= 0) return ndat;
    n = hasvalue(key) ? nemoinpr(getparam(key),dat,ndat) : 0;
    if (n < 0) error("Error %d parsing %s=%s",n,key,getparam(key));
    if (Qpixel) {
      for (i=0; i<n; i++) {
	dat[i] -= offset;
      }
    }

    return n;
}

void myprintf(string fmt,real v)
{
    printf(fmt,v);
    printf(" ");
}

int near_edge(int ix,int iy,int iz,int nx,int ny,int nz,int dim, int edge)
{
  if (ix<edge || ix > nx-1-edge) return 1;
  if (iy<edge || iy > ny-1-edge) return 1;
  if (dim==3)
    if (iz<edge || iz > nz-1-edge) return 1;
  return 0;
}
