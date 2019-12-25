/*
 * YAPP: Yet Another Plotting Package.
 *
 *	This implements all yapp calls as dummy stubs - they don't
 *	do anything, but just gets you going linking.
 *
 */

#include <stdinc.h>
#include <yapp.h>

int plinit(string pltdev, real xmin, real xmax, real ymin, real ymax)
{
 dprintf(0,"[YAPP_NULL: no graphics output]\n");
 return 0;
}

int plswap() 
{ return 0;}

real plxscale(real x, real y)
{ return 0;}

real plyscale(real x, real y)
{ return 0;}

int plltype(int lwid, int lpat)
{ return 0;}

int plline(real x, real y)
{ return 0;}

int plmove(real x, real y)
{ return 0;}

int plpoint(real x, real y)
{ return 0;}

int plcircle(real x, real y, real r)
{ return 0;}

int plcross(real x, real y, real s)
{ return 0;}

int plbox(real x, real y, real s)
{ return 0;}

int pljust(int jus)
{ return 0;}

int pltext(string msg, real x, real y, real hgt, real ang)
{ return 0;}

int plflush() 
{ return 0;}

int plframe()
{ return 0;}

int plstop()
{ return 0;}

void plcolor(int color)
{ }

int pl_matrix(real *frame,int nx,int ny,real xmin,real ymin,
	      real cell,real fmin,real fmax,real findex, real blankval)
{ return 0;}

int pl_contour(real *frame,int nx,int ny, int nc, real *c)
{ return 0;}

int pl_screendump(string fname)
{ return 0;}

int pl_getpoly(float *x,float *y,int n)
{ return 0;}

int pl_cursor(real *x,real *y, char *c)
{ return 0;}

int pl_interp(string cmd)
{ return 0;}
