/*
 * YAPP: Yet Another Plotting Package.
 *
 *	This implements all yapp calls as dummy stubs - they don't
 *	do anything, but just gets you going linking.
 */

#include <stdinc.h>

plinit(string pltdev, real xmin, real xmax, real ymin, real ymax)
{
 dprintf(0,"[YAPP_NULL: no graphics output]\n");
}

plswap() 
{}

real plxscale(real x, real y)
{ return 1.0;}

real plyscale(real x, real y)
{ return 1.0;}


plltype(int lwid, int lpat)
{}

plline(real x, real y)
{}

plmove(real x, real y)
{}

plpoint(real x, real y)
{}

plcircle(real x, real y, real r)
{}

plcross(real x, real y, real s)
{}

plbox(real x, real y, real s)
{}

pljust(int jus)
{}

pltext(string msg, real x, real y, real hgt, real ang)
{}

plflush() 
{}

plframe()
{}

plstop()
{}

plcolor(int color)
{}

pl_matrix(real *frame,int nx,int ny,real xmin,real ymin,
	  real cell,real fmin,real fmax,real findex)
{}

pl_contour(real *frame,int nx,int ny, int nc, real *c)
{}

pl_screendump(string fname)
{}

pl_getpoly(float *x,float *y,int n)
{
 return 0;
}

pl_cursor(float *x,float *y, char *c)
{
 return 0;
}

pl_interp(string cmd)
{}
