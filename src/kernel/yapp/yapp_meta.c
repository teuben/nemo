/*
 * YAPP: Yet Another Plotting Package.
 *
 *	This implements a simple ascii meta-file; can be edited
 *	Note that the units are in cm, and only accurate to 0.1mm
 *	Sometimes useful for debugging.
 */

#include <stdinc.h>

extern string yapp_string;

local stream yappstr = NULL;
local string yapp_dev;
local int ncolors = 16;      /* kind of faking pgplot */

int plinit(string pltdev, real xmin, real xmax, real ymin, real ymax)
{
    yappstr = stropen(yapp_string,"w");
    fprintf(yappstr,"plinit %5.2f %5.2f %5.2f %5.2f %s\n",
            xmin,xmax,ymin,ymax,yapp_string);
}

int plswap() 
{
    fprintf(yappstr,"plswap\n");
}

real plxscale(real x, real y)
{
    fprintf(yappstr,"plxscale %5.2f %5.2f\n",x,y);
}

real plyscale(real x, real y)
{
    fprintf(yappstr,"plyscale %5.2f %5.2f\n",x,y);
}

int plltype(int lwid, int lpat)
{
    fprintf(yappstr,"plltype %d %d\n",lwid,lpat);
}

int plline(real x, real y)
{
    fprintf(yappstr,"plline %5.2f %5.2f\n",x,y);
}

int plmove(real x, real y)
{
    fprintf(yappstr,"pldraw %5.2f %5.2f\n",x,y);
}

int plpoint(real x, real y)
{
    fprintf(yappstr,"plpoint %5.2f %5.2f\n",x,y);
}

int plcircle(real x, real y, real r)
{
    fprintf(yappstr,"plcircle %5.2f %5.2f %5.2f\n",x,y,r);    
}

int plcross(real x, real y, real s)
{
    fprintf(yappstr,"plcross %5.2f %5.2f %5.2f\n",x,y,s);
}

int plbox(real x, real y, real s)
{
    fprintf(yappstr,"plbox %5.2f %5.2f %5.2f\n",x,y,s);
}

int pljust(int jus)
{
    fprintf(yappstr,"pljust %d\n",jus);
}

int pltext(string msg, real x, real y, real hgt, real ang)
{
    fprintf(yappstr,"pltext %5.2f %5.2f %5.2f %5.2f %s\n",x,y,hgt,ang,msg);
}

int plflush() 
{
    fprintf(yappstr,"plflush\n");
}

int plframe()
{
    fprintf(yappstr,"plframe\n");
}

int plstop()
{
    fprintf(yappstr,"plstop\n");
    strclose(yappstr);
    yappstr=NULL;
}

int plncolors()
{
  return ncolors;
}

void plcolor(int color)
{
    fprintf(yappstr,"color %d\n",color);
}

void plpalette(real *r, real *g, real *b, int n)
{
    /* cannot handle them */
    fprintf(yappstr,"# pallete %d\n",n);
}


int pl_matrix(real *frame,int nx,int ny,real xmin,real ymin, real cell,
	  real fmin,real fmax,real findex, real blankval)
{
    fprintf(yappstr,"pl_matrix\n");
}

int pl_screendump(string fname)
{
    fprintf(yappstr,"pl_screendump\n");
}

int pl_getpoly(float *x, float *y, int n)
{
    fprintf(yappstr,"pl_getpoly\n");
    return 0;
}

int pl_interp(string cmd)
{
    fprintf(yappstr,"pl_interp\n");
}
