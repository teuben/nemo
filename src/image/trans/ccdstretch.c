/* 
 * CCDSTRETCH: stretch an image dimension
 *
 *       17-aug-2023     V0.1 written for Stuart's rosetta cube
 *       19-aug-2023     V0.5 added compress= and average=
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>
#include <string.h>

string defv[] = {
        "in=???\n       Input image cube",
	"out=???\n      Output image cube",
        "factor=2\n     Integer factor by which Z axis  is stretched",
	"axis=z\n       Axis to be stretched (z,y,z)",
	"compress=f\n   Compress (rebin) instead of stretch?",
	"average=t\n    Average when rebinning?",
	"VERSION=0.5\n  19-aug-2023 PJT",
	NULL,
};

string usage = "stretching (or compressing) an image cube along an axis";


#define CVI(x,y,z)  CubeValue(iptr,x,y,z)
#define CVO(x,y,z)  CubeValue(optr,x,y,z)


void nemo_main()
{
    stream   instr, outstr;
    imageptr iptr=NULL, optr;
    int      nx, ny, nz, nx2, ny2, nz2, i,j,k,k1,k2;
    int      axis = 0;  // valid are 1,2,3
    int      factor = getiparam("factor");
    real     rfactor, sum, aver;
    string   saxis = getparam("axis");
    bool     Qbin = getbparam("compress");
    bool     Qaver = getbparam("average");

    if (strchr("xX1", *saxis)) axis=1;
    if (strchr("yY2", *saxis)) axis=2;
    if (strchr("zZ3", *saxis)) axis=3;
    if (axis == 0) error("illegal axis %s (%d)", saxis, axis);

    rfactor = Qbin ?   1.0/factor   :  factor;
    aver = Qaver ?  factor : 1.0;
      

    instr = stropen(getparam("in"), "r");
    read_image( instr, &iptr);
    nx = Nx(iptr);	
    ny = Ny(iptr);
    nz = Nz(iptr);
    if (Qbin) {
      nx2 = axis==1 ? nx/factor : nx;
      ny2 = axis==2 ? ny/factor : ny;
      nz2 = axis==3 ? nz/factor : nz;
    } else {
      nx2 = axis==1 ? nx*factor : nx;
      ny2 = axis==2 ? ny*factor : ny;
      nz2 = axis==3 ? nz*factor : nz;
    }
    dprintf(0,"Input: %d x %d x %d   Output: %d x %d x %d   %s %s %s\n",
	    nx,ny,nz, nx2,ny2,nz2,
	    saxis, Qbin ? "compress" : "stretch ",  Qbin ?  (Qaver ? "aver" : "sum") : "");

    outstr = stropen(getparam("out"), "w");
    create_cube(&optr,nx2,ny2,nz2);
    copy_header(iptr,optr,1);
    /* fix cell size */
    if(axis==1) Dx(optr) /= rfactor;
    if(axis==2) Dy(optr) /= rfactor;
    if(axis==3) Dz(optr) /= rfactor;
    /* fix reference pixel */
    if(axis==1) Xref(optr) = (Xref(iptr)+0.5)*rfactor - 0.5;
    if(axis==2) Yref(optr) = (Yref(iptr)+0.5)*rfactor - 0.5;
    if(axis==3) Zref(optr) = (Zref(iptr)+0.5)*rfactor - 0.5;

    if (axis==3)                                 /* handle the Z axis */
      for (i=0; i<nx; i++) 
	for (j=0; j<ny; j++) 
	  if (Qbin) {
	    for (k=0, k1=0; k<nz2; k++)	{
	      for (k2=0, sum=0.0; k2<factor; k2++, k1++) {
		sum += CVI(i,j,k1);
	      }
	      CVO(i,j,k) = sum/aver;
	    }
	  } else
	    for (k=0, k1=0; k<nz; k++)	      
	      for (k2=0; k2<factor; k2++, k1++)
		CVO(i,j,k1) = CVI(i,j,k);
    else if (axis==2)                            /* handle the Y axis */
      for (i=0; i<nx; i++) 
	for (j=0; j<nz; j++) 
	  if (Qbin) {
	    for (k=0, k1=0; k<ny2; k++)	{
	      for (k2=0, sum=0.0; k2<factor; k2++, k1++) {
		sum += CVI(i,k1,j);
	      }
	      CVO(i,k,j) = sum/aver;
	    }
	  } else
	    for (k=0, k1=0; k<ny; k++)	      	    
	      for (k2=0; k2<factor; k2++, k1++)
		CVO(i,k1,j) = CVI(i,k,j);
    else if (axis==1)                            /* handle the X axis */
      for (i=0; i<ny; i++) 
	for (j=0; j<nz; j++) 
	  if (Qbin) {
	    for (k=0, k1=0; k<nx2; k++)	{
	      for (k2=0, sum=0.0; k2<factor; k2++, k1++) {
		sum += CVI(k1,i,j);
	      }
	      CVO(k,i,j) = sum/aver;
	    }
	  } else
	    for (k=0, k1=0; k<nx; k++)	      	    	    
	      for (k2=0; k2<factor; k2++, k1++)
		CVO(k1,i,j) = CVI(k,i,j);
    else
      error("bad axis %d", axis);
    write_image(outstr, optr);
}

