/* 
 * TWSPEED
 *
 *	30-mar-03  V1.0 resurrected from pspeed as twspeed
 *       3-aug-06  V1.1 optional h(Y) windowing (following TW84)
 *
OLD::
./apus1/teuben/nemo/usr/pjt/image/pspeed.c
./apus1/teuben/mars/teuben/d/image/pspeed.c
./apus1/nemo/usr/pjt/image/pspeed.c
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>
#include <image.h>
#include <moment.h> 

string defv[] = {
    "vel=???\n       Input mean velocity image",
    "den=???\n       Input density image",
    "xin=\n          Input mean position image",
    "pa=90\n         *** Position Angle of disk",
    "inc=0\n         *** Inclination of disk",
    "center=\n       *** Center of galaxy (map-center)",
    "step=1\n        Slit width in pixels",
    "window=0\n      Window functions (0=none, ..)",
    "VERSION=1.2\n   3-aug-06 PJT",
    NULL,
};

string usage="TW pattern speed of an object";

stream  vinstr, xinstr, instr;

imageptr viptr=NULL;			/* will be allocated dynamically */
imageptr iptr=NULL;			/* will be allocated dynamically */
imageptr xiptr=NULL;			/* will be allocated dynamically, if used */

int    nx,ny,nz,nsize;			/* actual size of map */
double xmin,ymin,zmin,dx,dy,dz;
double size;				/* size of frame (square) */
double cell;				/* cell or pixel size (square) */


void nemo_main(void)
{
  int  i, j, k, dir;
  real sum0, sum1, sum2, sum11,sum00, xmean, y;
    
  vinstr = stropen (getparam("vel"), "r");         /* velocity map */
  read_image (vinstr,&viptr);

  instr = stropen (getparam("den"), "r");         /* density map */
  read_image (instr,&iptr);

  if (hasvalue("xin")) {
    xinstr = stropen (getparam("xin"), "r");      /* optional Xmean map */
    read_image (xinstr,&xiptr);
  } 
  
  nx = Nx(iptr);	
  ny = Ny(iptr);
  xmin = Xmin(iptr);
  ymin = Ymin(iptr);
  dx = Dx(iptr);
  dy = Dy(iptr);

  dprintf(1,"# j sum0 sum1 sum2  sum1/sum0  sum2/sum0\n");
  for (j=0; j<ny; j++) {     /* loop over all lines parallel to the major axis */
    y = ymin + j*dy;
    sum0 = sum1 = sum2 = sum11 = sum00 = 0;
    for (i=0; i<nx; i++) {
      sum00 += MapValue(iptr,i,j);
      sum11 += (xmin + i*dx) * MapValue(iptr,i,j);
    }
    xmean = sum11/sum00;
    for (i=0; i<nx; i++) {          /* integrate TW quantities */
      sum0 += MapValue(iptr,i,j);                        /* integrate DEN map */
      sum1 += MapValue(viptr,i,j)*MapValue(iptr,i,j);    /* integrate VEL map */
      if (xiptr)
	sum2 += MapValue(xiptr,i,j); 
      else
	sum2 += xmean;
    }
    if (sum0 != 0)
      printf("%g %g  %g  %d %g %g %g     %g %g\n",
	     xmean, sum1/sum0, y,
	     j,sum0,sum1,sum2,sum1/sum0,sum2/sum0);

  }
  printf("# xmean vmean  y  j sum0 sum1 sum2  sum1/sum0  sum2/sum0\n");
}
