/* 
 * CCDZ: take a Z function of a cube
 *
 *	quick and dirty:  20-aug-2013
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
  "in=???\n       Input image file",
  "out=???\n      Output image file",
  "lut=$NEMODAT/rainbow1.lut\n     RGB response",
  "z=\n           z range, if used",
  "gamma=1\n      Gamma factor",
  "color=R,G,B\n  Pick R, G or B",
  "clip=\n        If used, clip values between -clip,clip or clip1,clip2",
  "VERSION=0.2\n  20-aug-2013 PJT",
  NULL,
};

string usage = "Z function of a cube";
string cvsid="$Id$";

local real peak_axis(imageptr iptr, int i, int j, int k, int axis);
local int  peak_find(int n, real *data, int *mask, int npeak);
local bool out_of_range(real);
local void image_oper(imageptr ip1, string oper, imageptr ip2);



void nemo_main()
{
    stream  instr, outstr;
    string  oper;
    int     i,j,k,nx, ny, nz, nx1, ny1, nz1;
    int     axis, mom;
    int     nclip, apeak, apeak1, cnt;
    imageptr iptr=NULL, iptr1=NULL, iptr2=NULL;      /* pointer to images */
    real    tmp0, tmp1, tmp2, tmp00, newvalue, peakvalue, scale, offset;
    real    *spec, ifactor, cv, clip[2];
    int     zrange[2];
    int     *smask;
    bool    Qclip = hasvalue("clip");

    instr = stropen(getparam("in"), "r");

    if (Qclip) {
      nclip = nemoinpr(getparam("clip"),clip,2);
      if (nclip<1) error("error parsing clip=%s",getparam("clip"));
      if (nclip==1) {
	clip[1] =  clip[0];
	clip[0] = -clip[1];
      }
    }

    nemoinpi(getparam("z"),zrange,2);
    dprintf(0,"Zrange: %d %d \n",zrange[0],zrange[1]);

    read_image( instr, &iptr);
    nx1 = nx = Nx(iptr);	
    ny1 = ny = Ny(iptr);
    nz1 = nz = Nz(iptr);

    nx1 = nx;   ny1 = ny;   nz1 = 1;
    spec = (real *) allocate(nz*sizeof(real));
    smask = (int *) allocate(nz*sizeof(int));
    dprintf(0,"Reducing %d*%d*%d to a %d*%d*%d cube\n",
	    nx,ny,nz, nx1,ny1,nz1);

    outstr = stropen(getparam("out"), "w");

    create_cube(&iptr1,nx1,ny1,nz1);
    create_cube(&iptr2,nx1,ny1,nz1);

    ifactor = 1.0;

    for(j=0; j<ny; j++)
      for(i=0; i<nx; i++) {
	tmp0 = 0.0;
	for(k=0; k<nz; k++) {
	  if(zrange[0]<=k && k<=zrange[1])
	    tmp0 += CubeValue(iptr,i,j,k);
	}
	CubeValue(iptr1,i,j,0) = tmp0;
      }
        

    Xmin(iptr1) = Xmin(iptr);
    Ymin(iptr1) = Ymin(iptr);
    Zmin(iptr1) = Zmin(iptr) + 0.5*(nz-1)*Dz(iptr);
    Dx(iptr1) = Dx(iptr);
    Dy(iptr1) = Dy(iptr);
    Dz(iptr1) = nz * Dz(iptr);
    
    Namex(iptr1) = Namex(iptr); /* care: we're passing a pointer */
    Namey(iptr1) = Namey(iptr);
    Namez(iptr1) = Namez(iptr);

    write_image(outstr, iptr1);
}


#define CV(p,i,j,k)   CubeValue(p,i,j,k)
#define LOOP(i,n)     for(i=0;i<n;i++)

local void image_oper(imageptr ip1, string oper, imageptr ip2)
{
  int i,j,k, nx,ny,nz;
  nx = Nx(ip1);
  ny = Ny(ip1);
  nz = Nz(ip1);
  if(nx!=Nx(ip2)) error("image_oper: NX size %d != %d",nx,Nx(ip2));
  if(ny!=Ny(ip2)) error("image_oper: NY size %d != %d",ny,Ny(ip2));
  if(nz!=Nz(ip2)) error("image_oper: NZ size %d != %d",nz,Nz(ip2));

  if (*oper== '+')
    LOOP(k,nz) LOOP(j,ny) LOOP(i,nx) CV(ip2,i,j,k) = CV(ip1,i,j,k) + CV(ip2,i,j,k);
  else if (*oper== '-')
    LOOP(k,nz) LOOP(j,ny) LOOP(i,nx) CV(ip2,i,j,k) = CV(ip1,i,j,k) - CV(ip2,i,j,k);
  else if (*oper== '*')
    LOOP(k,nz) LOOP(j,ny) LOOP(i,nx) CV(ip2,i,j,k) = CV(ip1,i,j,k) * CV(ip2,i,j,k);
  else if (*oper== '/')
    LOOP(k,nz) LOOP(j,ny) LOOP(i,nx) CV(ip2,i,j,k) = CV(ip1,i,j,k) / CV(ip2,i,j,k);
  else 
    error("invalid operator: %s",oper);
 
}
