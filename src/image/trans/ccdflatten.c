/* 
 * CCDFLATTEN: flattening an image of single pixel peaks, 
 *             a quick alternative to ccdmedian
 *
 *       2-aug-04       PJT     written
 *                      
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>
#include <moment.h>

string defv[] = {
        "in=???\n       Input image file",
	"out=???\n      Output image file",
	"nsigma=20\n	cutoff in terms of sigma",
	"n=1\n          width of area around pixel to get sigma",
	"VERSION=0.2\n  26-jan-2021 PJT",
	NULL,
};

string usage = "flattening an image of single pixel peaks";

string cvsid = "$Id$";


#if 0
#define CVI(x,y,z)  CubeValue(iptr,x,y,z)
#define CVO(x,y,z)  CubeValue(optr,x,y,z)
#else
#define CVI(x,y)    MapValue(iptr,x,y)
#define CVO(x,y)    MapValue(optr,x,y)
#endif

void nemo_main()
{
    stream  instr, outstr;
    int     nx, ny, nz, count, n;
    int     i,j,  i1, j1;
    real    ratio, nsigma, mean, sigma;
    imageptr iptr=NULL, optr;      /* pointer to images */
    Moment  m;

    nsigma = getdparam("nsigma");
    n = getiparam("n");

    instr = stropen(getparam("in"), "r");
    read_image( instr, &iptr);
    nx = Nx(iptr);	
    ny = Ny(iptr);
    nz = Nz(iptr);
    if (nz > 1) error("Cannot do 3D cubes properly; use 2D");

    outstr = stropen(getparam("out"), "w");
    create_cube(&optr,nx,ny,nz);
    Dx(optr) = Dx(iptr);
    Dy(optr) = Dy(iptr);
    Dz(optr) = Dz(iptr);
    Xmin(optr) = Xmin(iptr);
    Ymin(optr) = Ymin(iptr);
    Zmin(optr) = Zmin(iptr);
    Xref(optr) = Xref(iptr);
    Yref(optr) = Yref(iptr);
    Zref(optr) = Zref(iptr);
    Axis(optr) = Axis(iptr);

    count = 0;
    for (j=n; j<ny-n; j++) {
      for (i=n; i<nx-n; i++) {
	ini_moment(&m,2,0);
	for (j1=j-n; j1<=j+n; j1++) {
	  for (i1=i-n; i1<=i+n; i1++)
	    if (i1!=i || j1!=j)
	      accum_moment(&m, CVI(i1,j1), 1.0);
	}
	if (sigma_moment(&m) > 0) {
	  mean = mean_moment(&m);
	  sigma = sigma_moment(&m);
	  ratio = (CVI(i,j) - mean)/sigma;
	  ratio = ABS(ratio);
	  dprintf(1,"%d %d %g %g %g %d\n",i,j,CVO(i,j),mean,sigma,n_moment(&m));
	  if (ratio > nsigma) {
	    count++;
	    CVO(i,j) = mean_moment(&m);
	  } else
	    CVO(i,j) = CVI(i,j);
	} else
	  CVO(i,j) = CVI(i,j);
      }
    }
    dprintf(0,"ccdflatten: found %d/%d at nsigma=%g\n",count,nx*ny,nsigma);
    write_image(outstr, optr);
}

