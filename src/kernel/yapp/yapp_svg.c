/*
 * YAPP: Yet Another Plotting Package.
 *
 *      This implements a simple SVG (xml) plotting format
 *      See e.g. https://github.com/adishavit/simple-svg/blob/master/simple_svg_1.0.0.hpp
 *               https://github.com/davecom/SVGChart
 *
 *      Note that the plplot library has a pretty decent SVG output mode
 */

#include <stdinc.h>

extern string yapp_string;

local stream yappstr = NULL;
local string yapp_dev;
local int ncolors = 16;

local long svg_width  = 800;
local long svg_height = 800;

local real xc = 0;          /* the current pen position -- reminds me to COMPLOT */
local real yc = 0;

local real _xmin, _xmax, _ymin, _ymax;
local int _lineR=0, _lineG=0, _lineB=0;    /* colors are 0..255 */


local real x2w(real x) {
  return (x-_xmin)/(_xmax-_xmin)*svg_width;
}

local real y2h(real y) {
  return (y-_ymin)/(_ymax-_ymin)*svg_height;
}

int plinit(string pltdev, real xmin, real xmax, real ymin, real ymax)
{
  _xmin = xmin;
  _xmax = xmax;
  _ymin = ymin;
  _ymax = ymax;
    
  yappstr = stropen(yapp_string,"w");

  fprintf(yappstr,"<?xml version=\"1.0\" encoding=\"utf-8\"?>\n");
  fprintf(yappstr,"<svg version=\"1.1\" baseProfile=\"full\" width=\"%ld\" height=\"%ld\"",svg_width,svg_height);
  fprintf(yappstr," xmlns=\"http://www.w3.org/2000/svg\">\n");
}

int plswap() 
{
  //fprintf(yappstr,"plswap\n");
}

real plxscale(real x, real y)
{
  //fprintf(yappstr,"plxscale %5.2f %5.2f\n",x,y);
}

real plyscale(real x, real y)
{
  //fprintf(yappstr,"plyscale %5.2f %5.2f\n",x,y);
}

int plltype(int lwid, int lpat)
{
  //fprintf(yappstr,"plltype %d %d\n",lwid,lpat);
}

int plline(real x, real y)
{
    fprintf(yappstr,"<line x1=\"%g\" y1=\"%g\"  x2=\"%g\" y2=\"%g\"",
	    x2w(xc),y2h(yc),x2w(x),y2h(y));
    fprintf(yappstr," stroke=\"rgb(%d,%d,%d)\"", _lineR, _lineG, _lineB);
    fprintf(yappstr,"/>\n");
    xc = x;
    yc = y;
}

int plmove(real x, real y)
{
  xc = x;
  yc = y;
  //fprintf(yappstr,"pldraw %5.2f %5.2f\n",x,y);
}

int plpoint(real x, real y)
{
  //fprintf(yappstr,"plpoint %5.2f %5.2f\n",x,y);
}

int plcircle(real x, real y, real r)
{
  //fprintf(yappstr,"plcircle %5.2f %5.2f %5.2f\n",x,y,r);    
}

int plcross(real x, real y, real s)
{
  //fprintf(yappstr,"plcross %5.2f %5.2f %5.2f\n",x,y,s);
}

int plbox(real x, real y, real s)
{
  //fprintf(yappstr,"plbox %5.2f %5.2f %5.2f\n",x,y,s);
}

int pljust(int jus)
{
  //fprintf(yappstr,"pljust %d\n",jus);
}

#if 0
void SVGPainter::DrawText (int inX, int inY, const char *inString) {
        svgContent << "<text x=\"" << inX << "\" y=\"" << inY;
        svgContent << "\" fill=\"" << "rgb(" << lineRed << ",";
        svgContent << lineGreen << "," << lineBlue << ")" << "\">\n";
        svgContent << inString << "</text>\n";
    }

void SVGPainter::DrawRotatedText (int inX, int inY, float inDegrees, const char *inString) {
        svgContent << "<text x=\"" << inX << "\" y=\"" << inY << "\" transform=";
        svgContent << "\"rotate(" << inDegrees;
        svgContent << "," << inX << "," << inY << ")\">";
        svgContent << inString << "</text>\n";
    }
#endif

int pltext(string msg, real x, real y, real hgt, real ang)
{
  fprintf(yappstr,"<text x=\"%g\" y=\"%g\" transform=\"rotate(%g,%g,%g)\"> %s</text>\n",
	  x2w(x),y2h(y), ang, x2w(x),y2h(y),msg);
}

int plflush() 
{
  //    fprintf(yappstr,"plflush\n");
}

int plframe()
{
  //fprintf(yappstr,"plframe\n");
}

int plstop()
{
    fprintf(yappstr,"</svg>\n");
    strclose(yappstr);
    yappstr=NULL;
}

int plncolors()
{
  return ncolors;
}

void plcolor(int color)
{
  //fprintf(yappstr,"color %d\n",color);
}

void plpalette(real *r, real *g, real *b, int n)
{
    /* cannot handle them */
    //fprintf(yappstr,"# pallete %d\n",n);
}


int pl_matrix(real *frame,int nx,int ny,real xmin,real ymin, real cell,
	  real fmin,real fmax,real findex, real blankval)
{
  //fprintf(yappstr,"pl_matrix\n");
}

int pl_screendump(string fname)
{
  //fprintf(yappstr,"pl_screendump\n");
}

int pl_getpoly(float *x, float *y, int n)
{
  //fprintf(yappstr,"pl_getpoly\n");
    return 0;
}

int pl_interp(string cmd)
{
  //fprintf(yappstr,"pl_interp\n");
}
