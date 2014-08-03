/* 
 * CCDPATCH: patch a value in a polygon region
 *
 *       2-aug-2014     quick hack for DataBay_MD               pjt
 *
 *  @todo   allow circles, boxes, etc.
 *                      
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
  "in=???\n       Input image file",
  "out=???\n      Output image file",
  "patch=0.0\n    Data value to patch in",
  "polygon=\n     Polygon",
  "box=\n         Box xmin,ymin,xmax,ymax ** not implemented **",
  "circle=\n      XCenter,Ycenter,Radius ** not implemented **",
  "ellipse=\n     XCenter,Ycenter,Major,Minor,PositionAngle ** not implemented **",
  "VERSION=0.2\n  2-aug-2014 PJT",
  NULL,
};

string usage = "patch up regions from polygon/shapes";

string cvsid = "$Id$";


#define MAXPOLY       128

local int    np;
local real   polygon[2*MAXPOLY];
local real   xp[MAXPOLY], yp[MAXPOLY];

int  pairs(string p, real *pairs, real *x, real *y, int maxp);
int  inpolygon (int n, real *x, real *y, real x0, real y0);


void nemo_main(void)
{
  stream   instr, outstr;
  int      m, n, nx, ny, nz;             /* size of map */
  int      i,j,k, di, dj, npatch;
  imageptr iptr=NULL;                    /* pointer to images */
  real     patch = getdparam("patch");
  real     xi, yj;

  np = pairs(getparam("polygon"),polygon,xp,yp,2*MAXPOLY);

  instr = stropen(getparam("in"), "r");
  
  read_image( instr, &iptr);
  
  nx = Nx(iptr);	
  ny = Ny(iptr);
  nz = Nz(iptr);
  
  outstr = stropen(getparam("out"), "w");
  
  npatch = 0;
  for (k=0; k<nz; k++) {
    for (j=0; j<ny; j++) {
      yj = j+1;
      for (i=0; i<nx; i++) {
	xi = i+1;
	if (inpolygon(np,xp,yp,xi,yj)) {
	  npatch++;
	  CubeValue(iptr,i,j,k) = patch;
	}
      }
    }
    dprintf(0,"Found %d values to patch\n",npatch);
  } /* k */
  write_image(outstr, iptr);
}

/* the code below copied from tabpolygon, unmodified */

int pairs(string p, real *pairs, real *xp, real *yp, int maxp)
{
  int n, i;

  n = nemoinpr(p,pairs,maxp);
  if (n>=0 && n<6) error("%s needs at least 6 values for polygon=",p);
  if (n%2) error("Need even number of points for x,y pairs in polygon=");
  if (n<0) error("Parsing error %s",p);
  for (i=0; i<n; i++) {
    if (i%2)
      yp[i/2] = pairs[i];
    else
      xp[i/2] = pairs[i];
  }
  n = n/2;
  for (i=0; i<n; i++) 
    dprintf(1,"polygon-%d: %g %g\n",i,xp[i],yp[i]);
  return n;
}


/*	Early version, for updates see snapplot.c or some library */
/*  The current algorithm has some flaws, may not be the fastest,
 *  (I have seen a faster one somewhere, but can't recall where),
 *  works. It's from CACM aug 1962, algorithm 112 by M. Shimrat (p434)
 *  but see also remark in CACM dec 1962, p606)
 */

/*		Code for this is now in editplot/snapplot */

int inpolygon_old (int n, real *x, real *y, real x0, real y0)
{
  int i,b;              /* b=0:outside b=1:inside */
        
  x[n] = x[0];        /* make sure polygon is closed */
  y[n] = y[0];

  b = 0;              /* set initially to false */
  for (i=0; i<n; i++) {
    if ( ((y0<y[i])==(y0>y[i+1])) &&
	 ((x0-x[i]-(y0-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i])) < 0 ))
      b = !b;
  }
  return b;
}

int inpolygon(int nvert, real *vertx, real *verty, real testx, real testy)
{
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
     (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }
  return c;
}

