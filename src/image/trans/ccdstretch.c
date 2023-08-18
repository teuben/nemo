/* 
 * CCDSTRETCH: stretch image dimensions via interpolation
 *             opposite of CCDSUB
 *
 *       17-aug-2023       PJT     written for Stuart's retirement cube
 *                      
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
        "in=???\n       Input image cube",
	"out=???\n      Output image cube",
        "factor=2\n     Integer factor by which Z axis  is stretched",
	"axis=z\n       Axis to be stretched (z,y,z)",
	"VERSION=0.2\n  18-aug-2023 PJT",
	NULL,
};

string usage = "stretching an image cube along an axis";



#define CVI(x,y,z)  CubeValue(iptr,x,y,z)
#define CVO(x,y,z)  CubeValue(optr,x,y,z)

void nemo_main()
{
    stream  instr, outstr;
    int     nx, ny, nz, nx2, ny2, nz2;
    int     i,j, k, k1, k2;
    int     axis = 0;
    int     factor;
    imageptr iptr=NULL, optr;      /* pointer to images */
    string   saxis = getparam("axis");

    factor = getiparam("factor");
    if (*saxis == 'x') axis=1;
    if (*saxis == 'y') axis=2;
    if (*saxis == 'z') axis=3;
    if (axis == 0) error("illegal axis %d", axis);
    dprintf(0,"axis=%s  %d\n",saxis,axis);

    instr = stropen(getparam("in"), "r");
    read_image( instr, &iptr);
    nx = Nx(iptr);	
    ny = Ny(iptr);
    nz = Nz(iptr);
    nx2 = axis==1 ? nx*factor : nx;
    ny2 = axis==2 ? ny*factor : ny;
    nz2 = axis==3 ? nz*factor : nz;
    dprintf(0,"Input: %d x %d x %d   Output: %d x %d x %d \n",
	    nx,ny,nz, nx2,ny2,nz2);

    outstr = stropen(getparam("out"), "w");
    create_cube(&optr,nx2,ny2,nz2);
    copy_header(iptr,optr,1);
    /* fix cell size */
    if(axis==1) Dx(optr) /= factor;
    if(axis==2) Dy(optr) /= factor;
    if(axis==3) Dz(optr) /= factor;
    /* fix reference pixel, recalling integer pixel reference is in center of pixel */
    if(axis==3) Xref(optr) = (Xref(iptr)+0.5)*factor - 0.5;
    if(axis==3) Yref(optr) = (Yref(iptr)+0.5)*factor - 0.5;
    if(axis==3) Zref(optr) = (Zref(iptr)+0.5)*factor - 0.5;

    if (axis==3) {
      for (j=0; j<ny; j++) {
	for (i=0; i<nx; i++) {
	  for (k=0, k1=0; k<nz; k++) {
	    for (k2=0; k2<factor; k2++, k1++) {
	      CVO(i,j,k1) = CVI(i,j,k);
	    }
	  }
	}
      }
    } else if (axis==2) {
      for (j=0; j<nz; j++) {
	for (i=0; i<nx; i++) {
	  for (k=0, k1=0; k<ny; k++) {
	    for (k2=0; k2<factor; k2++, k1++) {
	      CVO(i,k1,j) = CVI(i,k,j);
	    }
	  }
	}
      }
    } else if (axis==1) {
      for (j=0; j<nz; j++) {
	for (i=0; i<ny; i++) {
	  for (k=0, k1=0; k<nx; k++) {
	    for (k2=0; k2<factor; k2++, k1++) {
	      CVO(k1,i,j) = CVI(k,i,j);
	    }
	  }
	}
      }
    } else
      error("bad axis %d", axis);
    MapMin(optr) = MapMin(iptr);
    MapMax(optr) = MapMax(iptr);
    dprintf(0,"New data min/max: %g %g\n",MapMin(optr),MapMax(optr));
    write_image(outstr, optr);
}

