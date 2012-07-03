/* 
 * CCDRT: apply some radiative transfer along the Z axis to create an XY map
 *
 *	quick and dirty, derived from ccdmom:  3-jul-2012

Some reading material:

http://www.cv.nrao.edu/course/astr534/Radxfer.html
http://en.wikipedia.org/wiki/Radiative_transfer

 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
        "in=???\n       Input intensity image file",
	"out=???\n      Output image file",
	"peak=f\n       Use peak value",
	"VERSION=0.1\n  3-jul-2012 PJT",
	NULL,
};

string usage = "apply some radiative transfer in a cube to create a map";
string cvsid="$Id$";

local real peak_axis(imageptr iptr, int i, int j, int k);
local bool out_of_range(real);



void nemo_main()
{
    stream  instr, outstr;
    int     nx, ny, nz, nx1, ny1, nz1;
    int     axis, mom;
    int     i,j,k, apeak, cnt;
    imageptr iptr=NULL, iptr1=NULL, iptr2=NULL;      /* pointer to images */
    real    tmp0, tmp1, tmp2, tmp00, newvalue, peakvalue, scale, offset;
    bool    Qpeak;

    instr = stropen(getparam("in"), "r");
    mom = 0;
    axis = 3;
    Qpeak = getbparam("peak");

    read_image( instr, &iptr);

    nx1 = nx = Nx(iptr);	
    ny1 = ny = Ny(iptr);
    nz1 = 1;
    nz  = Nz(iptr);

    outstr = stropen(getparam("out"), "w");

    create_cube(&iptr1,nx1,ny1,nz1);
    create_cube(&iptr2,nx1,ny1,nz1);


    scale = Dz(iptr);
    offset = Zmin(iptr);
    for(j=0; j<ny; j++)
      for(i=0; i<nx; i++) {
	tmp0 = tmp00 = tmp1 = tmp2 = 0.0;
	cnt = 0;
	peakvalue = CubeValue(iptr,i,j,0);
	for(k=0; k<nz; k++) {
	  if (out_of_range(CubeValue(iptr,i,j,k))) continue;
	  cnt++;
	  tmp0 += CubeValue(iptr,i,j,k);
	  tmp00 += sqr(CubeValue(iptr,i,j,k));
	  if (CubeValue(iptr,i,j,k) > peakvalue) {
	    apeak = k;
	    peakvalue = CubeValue(iptr,i,j,k);
	  }
	}
	if (cnt==0 || tmp0==0.0) {
	  newvalue = 0.0;
	} else {
	  if (Qpeak) 
	    newvalue = peakvalue;
	  else
	    newvalue = tmp0;
	}
	CubeValue(iptr1,i,j,1) = newvalue;
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


/*
 * return location of peak for (-1,y1) (0,y2) (1,y3)
 * as determined from the 2d polynomial going through
 * these 3 points
 */ 

local real peak_axis(imageptr iptr, int i, int j, int k)
{
  real y1, y2, y3;
  y2 = CubeValue(iptr,i,j,k);
  k -= 1;
  y1 = CubeValue(iptr,i,j,k);
  k += 2;
  y3 = CubeValue(iptr,i,j,k);

  if (y1+y3 == 2*y2) return 0.0;
  return 0.5*(y1-y3)/(y1+y3-2*y2);
}


local bool out_of_range(real x)
{
  return FALSE;
}
