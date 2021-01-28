/* 
 * CCDGAUSSFIT: fit a gaussian to a point and remove it
 *
 *      30-jul-04       PJT     written
 *                      
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>
#include <lsq.h>

string defv[] = {
  "in=???\n       Input image file",
  "out=???\n      Output image file",
  "x=\n           Initial guess of X position (0..nx-1)",
  "y=\n           Initial guess of Y position (0..ny-1)",
  "r1=\n          Inner radius (in pixels) of gaussian",
  "r2=\n          Outer radius (in pixels) for background estimation",
  "n=2\n          Order of 2D polynomial for background fit",
  "VERSION=0.2\n  26-jan-2021 PJT",
  NULL,
};

string usage = "fit a gaussian to a point and remove it";
  


#if 0
#define CVI(x,y,z)  CubeValue(iptr,x,y,z)
#define CVO(x,y,z)  CubeValue(optr,x,y,z)
#else
#define CVI(x,y)    MapValue(iptr,x,y)
#define CVO(x,y)    MapValue(optr,x,y)
#endif


void nemo_main(void)
{
  stream  instr, outstr;
  int     nx, ny, nz, mode, badval;
  int     i,j,k, n, n1, i0,i1, j0,j1;
  int     ix, iy, count;
  real    r, r1, r2, mapval, gauval, dx, dy;
  real    mat1[3*3], vec1[3], sol1[3], a1[4];
  real    mat2[5*5], vec2[5], sol2[5], a2[6];
  imageptr iptr=NULL, optr;      /* pointer to images */
  
  instr = stropen(getparam("in"), "r");
  read_image( instr, &iptr);
  nx = Nx(iptr);	
  ny = Ny(iptr);
  nz = Nz(iptr);
  if (nz > 1) error("Cannot do 3D cubes properly; use 2D");
  
  ix = getiparam("x");
  iy = getiparam("y");
  r1 = getdparam("r1");
  r2 = getdparam("r2");

  
  outstr = stropen(getparam("out"), "w");
  create_cube(&optr,nx,ny,nz);
  Dx(optr) = Dx(iptr);
  Dy(optr) = Dy(iptr);
  Dz(optr) = Dz(iptr);
  Xmin(optr) = Xmin(iptr);
  Ymin(optr) = Ymin(iptr);
  Zmin(optr) = Zmin(iptr);

  i0 = ix-r2-1;                 /* i0..i1 is the X range around star */
  if (i0<0) i0=0;
  i1 = ix+r2+1;
  if (i1>=nx) i1=nx-1;

  j0 = iy-r2-1;                 /* j0..j1 is the Y range around star */
  if (j0<0) j0=0;
  j1 = iy+r2+1;
  if (j1>=ny) j1=ny-1;

  dprintf(0,"Xrange: %d - %d    Yrange: %d - %d\n",i0,i1,j0,j1);
  dprintf(0,"r1=%g r2=%g\n",r1,r2);

  /*  loop over image accumulating data for background fit between r1 and r2 */
  if (r1 < r2) {
    lsq_zero(3,mat1,vec1);
    count = 0;
    for (j=j0; j<j1; j++) {
      dy = j-iy;
      for (i=i0; i<i1; i++) {
	dx = i-ix;
	r = sqrt(dx*dx+dy*dy);
	if (r<r1 || r>r2) continue;
	/* accumulate: */
	dprintf(1,"plane: x,y=%d %d r=%g\n",i,j,r);
	count++;
	a1[0] = 1.0;
	a1[1] = dx;
	a1[2] = dy;
	a1[3] = MapValue(iptr,i,j);
	lsq_accum(3,mat1,vec1,a1,1.0);
      }
    }
    dprintf(0,"%d points for plane fit\n",count);
    lsq_solve(3,mat1,vec1,sol1);
  } else {
    dprintf(0,"No planar fit\n",count);
    sol1[0] = sol1[1] = sol1[2] = 0;
  }


  /*  loop over image subtracting background between 0 and r2 */
  /*  also accumulate data for the gaussfit between 0 and r1 */
  lsq_zero(5,mat2,vec2);
  badval=0;
  count=0;
  for (j=j0; j<j1; j++) {
    dy = j-iy;
    for (i=i0; i<i1; i++) {
      dx = i-ix;
      r = sqrt(dx*dx+dy*dy);
      if (r>r1) continue;
      /* accumulate the gaussian data, but subtract the planar fit first */
      /* discard if data negative, user perhaps choose r1 too large      */
      count++;
#if 1
      mapval = MapValue(iptr,i,j) - sol1[0] - sol1[1]*dx - sol1[2]*dy;
#else
      mapval = MapValue(iptr,i,j);
#endif
      if (mapval <= 0.0) {
	dprintf(1,"mapval < 0::  x,y=%d %d mapval=%g %g\n",i,j,mapval,MapValue(iptr,i,j));
	badval++;
	continue;
      }
      dprintf(1,"gauss:  x,y=%d %d r=%g mapval=%g\n",i,j,r,mapval);
      a2[0] = 1.0;
      a2[1] = dx;
      a2[2] = dy;
      a2[3] = dx*dx;
      a2[4] = dy*dy;
      a2[5] = log(mapval);
      lsq_accum(5,mat2,vec2,a2,1.0);
    }
  }
  lsq_solve(5,mat2,vec2,sol2);
  dprintf(0,"%d points for gauss fit\n",count);
  dprintf(0,"Fit: %g %g %g %g %g \n",sol2[0],sol2[1],sol2[2],sol2[3],sol2[4]);
  if (1) {
    real xcen,ycen,sigx,sigy,amp;
    if (sol2[3] > 0 || sol2[4] > 0) error("bad gaussfit");
    sigx = sqrt(-0.5/sol2[3]);
    sigy = sqrt(-0.5/sol2[4]);
    xcen = sol2[1]/(sigx*sigx);
    ycen = sol2[2]/(sigy*sigy);
    amp = exp(sol2[0] + xcen*xcen/(2*sigx*sigx) + ycen*ycen/(2*sigy*sigy));
    dprintf(0,"Fit: a=%g x0=%g y0=%g sigx=%g sigy=%g \n",amp,xcen,ycen,sigx,sigy);
    
  }
  if (badval)
    warning("%d points around (%d,%d,%g,%g) have negative intensity after plane subtraction",
	    badval,ix,iy,r1,r2);

  /*  loop over image subtracting gaussfit between 0 and r1 */
  for (j=j0; j<j1; j++) {
    dy = j-iy;
    for (i=i0; i<i1; i++) {
      dx = i-ix;
      r = sqrt(dx*dx+dy*dy);
      if (r<r1) {
	gauval = sol2[0] + sol2[1]*dx + sol2[2]*dy + sol2[3]*dx*dx + sol2[4]*dy*dy;
#if 1
	MapValue(optr,i,j) = MapValue(iptr,i,j) - exp(gauval);
#else
	MapValue(optr,i,j) = exp(gauval);
#endif
      } else
	MapValue(optr,i,j) = MapValue(iptr,i,j);
    }
  }

  /* write out resulting image */
  write_image(outstr, optr);
}



/*
 *  both of these need a fitting routine of the kind (series cut by n=)
 *     f(x,y) = f0 +                                                        n=0
 *              f1.x +f2.y +                                                n=1
 *              f11.x.x + f12.x.y + f22.y.y +                               n=2
 *              f111.x.x.x + f112.x.x.y + f122.x.y.y + f222.y.y.y + ....    n=3
 *
 */
