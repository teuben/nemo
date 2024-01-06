/*
 * YAPP: Yet Another Plotting Package.
 *
 *      This uses the C interface to GNUPLOT using unix pipes (originally
 *      written by Nicolas Devillard) - 
 *      See https://github.com/longradix/gnuplot_i for the most recent
 *      implementation
 *      
 *      NOTE: this version is the placeholder until actually implemented
 */

#include <stdinc.h>
#include "gnuplot_i.c"       // needs to be installed in $NEMOLIB

extern string yapp_string;

//local stream yappstr = NULL;
//local string yapp_dev;
local real xp, yp;
local int ncolors = 16;      /* kind of faking pgplot */

local gnuplot_ctrl *handle;
local char cmd[128];

int plinit(string pltdev, real xmin, real xmax, real ymin, real ymax)
{
  dprintf(0,"yapp_gnuplot\n");
  handle = gnuplot_init();
  if (1) {
    gnuplot_setterm(handle, "wxt", 600, 600);  // 'enhanced font 'Verdana,10' persist'
  } else {
    gnuplot_cmd(handle, "set terminal png size 800, 800");
    gnuplot_cmd(handle, "set output \"yapp.png\"");
  }
  //gnuplot_setterm(handle, "dumb", 150, 40);
  
  //yappstr = stropen(yapp_string,"w");
  //fprintf(yappstr,"plinit %5.2f %5.2f %5.2f %5.2f %s\n",xmin,xmax,ymin,ymax,yapp_string);
  gnuplot_cmd(handle, "set xrange[0:20]");
  gnuplot_cmd(handle, "set yrange[0:20]");
  // Styling 
  gnuplot_cmd(handle, "set border linewidth 1.5");
  gnuplot_cmd(handle, "set pointsize 1.5");
  //  1=+ 2=x   3=x+  4=open sq
  gnuplot_cmd(handle, "set style line 1 lc rgb '#0060ad' pt 5");   // square
  gnuplot_cmd(handle, "set style line 2 lc rgb '#0060ad' pt 7");   // circle
  gnuplot_cmd(handle, "set style line 3 lc rgb '#0060ad' pt 9");   // triangle 
  xp = yp = 0.0;

  print_gnuplot_handle(handle);
  return 0;
}

int plswap() 
{
  return 0;
}

real plxscale(real x, real y)
{
  return 0;
}

real plyscale(real x, real y)
{
  return 0;
}

int plltype(int lwid, int lpat)
{
  return 0;
}

int plline(real x, real y)
{
  return 0;
}

int plmove(real x, real y)
{
  xp = x;
  yp = y;
  return 0;
}

int plpoint(real x, real y)
{
  gnuplot_cmd(handle, "plot '-' w p ls 1");
  sprintf(cmd, "%f %f", x, y);
  gnuplot_cmd(handle, cmd);
  gnuplot_cmd(handle, "e");

  dprintf(0,"plpoint: %s\n",cmd);
	     
  //gnuplot_cmd(handle, "set \n");
  return 0;
}

int plcircle(real x, real y, real r)
{
  return 0;
}

int plcross(real x, real y, real s)
{
  return 0;
}

int plbox(real x, real y, real s)
{
  return 0;
}

int pljust(int jus)
{
  return 0;
}

int pltext(string msg, real x, real y, real hgt, real ang)
{
  return 0;
}

int plflush() 
{
  return 0;
}

int plframe()
{
  gnuplot_resetplot(handle);
  return 0;
}

int plstop()
{
  printf("Press Enter to continue\n");
  while (getchar() != '\n') {}
  gnuplot_close(handle);
  return 0;
}

int plncolors()
{
  return ncolors;
}

void plcolor(int color)
{

}

void plpalette(real *r, real *g, real *b, int n)
{
}


int pl_matrix(real *frame,int nx,int ny,real xmin,real ymin, real cell,
	  real fmin,real fmax,real findex, real blankval)
{
  return 0;
}

int pl_screendump(string fname)
{
  return 0;
}

int pl_getpoly(float *x, float *y, int n)
{
  return 0;
}

int pl_interp(string cmd)
{
  return 0;
}

int pl_cursor(float *x, float *y, char *c)
{
  warning("no pl_cursor");
  return 0;
}
