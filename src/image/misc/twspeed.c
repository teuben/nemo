/* 
 * TWSPEED
 *
 *	30-mar-03  V1.0 resurrected from pspeed as twspeed
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
    "center=\n       *** Center of galaxy (map-center)",
    "step=1\n        Slit width in pixels",
    "VERSION=1.1\n   22-apr-04 PJT",
    NULL,
};

string usage="TW pattern speed of an object";

stream  vinstr, xinstr, instr;

imageptr viptr=NULL;			/* will be allocated dynamically */
imageptr xiptr=NULL;			/* will be allocated dynamically */
imageptr iptr=NULL;			/* will be allocated dynamically */


int    nx,ny,nz,nsize;			/* actual size of map */
double xmin,ymin,zmin,dx,dy,dz;
double size;				/* size of frame (square) */
double cell;				/* cell or pixel size (square) */



void nemo_main(void)
{
  int  i, j, k, dir;
  real sum0, sum1, sum2, sum11,sum00, xmean;
    
  vinstr = stropen (getparam("vel"), "r");         /* velocity map */
  read_image (vinstr,&viptr);
  strclose(vinstr);

  if (hasvalue("xin")) {
    xinstr = stropen (getparam("xin"), "r");      /* optional Xmean map */
    read_image (xinstr,&xiptr);
    strclose(xinstr);
  } else
    xinstr = NULL;
  
  instr = stropen (getparam("den"), "r");         /* density map */
  read_image (instr,&iptr);
  strclose(instr);

  
  nx = Nx(iptr);	
  ny = Ny(iptr);
  xmin = Xmin(iptr);
  ymin = Ymin(iptr);
  dx = Dx(iptr);
  dy = Dy(iptr);
  

  printf("# j sum0 sum1 sum2  sum1/sum0  sum2/sum0\n");
  for (j=0; j<ny; j++) {     /* loop over all lines parallel to the major axis */
    sum0 = sum1 = sum2 = sum11 = sum00 = 0;
    for (i=0; i<nx; i++) {
      sum00 += MapValue(iptr,i,j);
      sum11 += (xmin + i*dx) * MapValue(iptr,i,j);
    }
    xmean = sum11/sum00;
    for (i=0; i<nx; i++) {          /* integrate TW quantaties */
      sum0 += MapValue(iptr,i,j);                        /* integrate DEN map */
      sum1 += MapValue(viptr,i,j)*MapValue(iptr,i,j);    /* integrate VEL map */
      if (xiptr)
	sum2 += MapValue(xiptr,i,j); 
      else
	sum2 += xmean;
    }
    if (sum0 != 0)
      printf("%g %g   %d %g %g %g     %g %g\n",
	     xmean, sum1/sum0,
	     j,sum0,sum1,sum2,sum1/sum0,sum2/sum0);

  }
  printf("# xmean vmean   j sum0 sum1 sum2  sum1/sum0  sum2/sum0\n");
}
