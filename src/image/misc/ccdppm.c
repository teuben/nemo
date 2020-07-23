/* 
 * CCDPPM:   create a PPM from a CCD image
 *
 *	17-dec-03  V1.0 created   (Wright's 100 year first flight anniversary!)
 *                  1.1 added 8bit=
 *
 *  TODO
 *    - handle LUT's that are not 256 in size (e.g. ds9 uses 200 mostly)
 *    - handle lut='sls' to '$NEMODAT/lut/sls.lut' type translation
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>
#include <image.h>
#include <table.h>

string defv[] = {
    "in=???\n       Input image filename",
    "out=???\n      Output PPM filename",
    "min=\n         Minimum override",
    "max=\n         Maximum override",
    "bad=\n         Use this as masking value to ignore data",
    "power=1\n      Gamma factor applied to input data",
    "mean=\n        Mean value for ASINH scaling",
    "sigma=\n       Sigma value for ASINH scaling",
    "lut=\n         optional selection of a lut from NEMODAT/lut",
    "8bit=t\n       8bit?  (or 24bit) [not implemented]",
    "VERSION=1.3\n  19-jul-2020 PJT",
    NULL,
};

string usage="convert image to a PPM file, with optional color cable";

string	infile, outfile;       		/* file names */
stream  instr, outstr;			/* file streams */
imageptr iptr=NULL;			/* will be allocated dynamically */
int    nx,ny,nz;			/* actual size of map */

#define MAXCOL   256

int  ncol;
real  Red[MAXCOL],  Green[MAXCOL],  Blue[MAXCOL];
real dRed[MAXCOL], dGreen[MAXCOL], dBlue[MAXCOL];

void get_rgb(real data, real min, real max, int *R, int *G, int *B, bool);
int  get_lut(string, int, real*, real*, real*);

void nemo_main()
{
  int  i, j, red, green, blue;
  real x, xmin, xmax, bad, low, power, mean, sigma, slope;
  bool Qmin, Qmax, Qbad, Q8bit, Qasinh;
  
  instr = stropen (getparam("in"), "r");
  read_image (instr,&iptr);
  strclose(instr);
  outstr = stropen(getparam("out"), "w");
  Q8bit = getbparam("8bit");
  power = getrparam("power");
  Qasinh = hasvalue("mean") && hasvalue("sigma");
  if (Qasinh) {
    mean = getrparam("mean");
    sigma = getrparam("sigma");
  }

  if (!Q8bit) warning("24bit not working yet");
  
  nx = Nx(iptr);	
  ny = Ny(iptr);
  Qmin = hasvalue("min");
  if (Qmin) xmin = getdparam("min");
  Qmax = hasvalue("max");
  if (Qmax) xmax = getdparam("max");
  Qbad = hasvalue("bad");
  if (Qbad) bad = getdparam("bad");

  if (hasvalue("lut"))
    ncol = get_lut(getparam("lut"),MAXCOL,Red,Green,Blue);
  else {
    ncol = 256;
    for (i=0; i<ncol; i++) {
      Red[i] = Green[i] = Blue[i] = i/255.0;
      dRed[i] = dGreen[i] = dBlue[i] = 1.0/255.0;
    }
  }

  if (!Qmin || !Qmax) {                      /* loop over data to get xmin, xmax */
    if (!Qmin) xmin = MapValue(iptr,0,0);
    if (!Qmax) xmax = MapValue(iptr,0,0);
    for (j=0; j<ny; j++) {
      for (i=0; i<nx; i++) {
	x = MapValue(iptr,i,j);
	if (!Qmin && x<xmin) xmin=x;
	if (!Qmax && x>xmax) xmax=x;
      }
    }
    dprintf(0,"Auto-scaling, xmin=%g xmax=%g\n",xmin,xmax);
  }
  
  if (Qbad) {                              /* loop over data to flag bad data */
    dprintf(0,"Flag bad data\n");
    if (xmin<0)
      low = 2*xmin;
    else if (xmin>0)
      low = 0.5*xmin;
    else
      low = -xmax;
    for (j=0; j<ny; j++) {
      for (i=0; i<nx; i++) {
	x =  MapValue(iptr,i,j);
	if (Qbad && x==bad) MapValue(iptr,i,j) = low;
      }
    }
  }

  if (Qasinh) {                          /* if requested, apply the ASINH scaling */
    dprintf(0,"Using ASINH scaling with mean=%g sigma=%g\n",mean,sigma);
    slope = 1.0/asinh(xmax-xmin);
    for (j=0; j<ny; j++) {
      for (i=0; i<nx; i++) {
	x =  MapValue(iptr,i,j);
	MapValue(iptr,i,j) = asinh((x-mean)/sigma);
      }
    }
    xmin = 0;
    xmax = 1;
  }

  if (power != 1.0) {                   /* if requested, apply a gamma factor */
    dprintf(0,"Using gamma factor power=%g scaling\n",power);
    slope = 1.0/(xmax-xmin);    
    for (j=0; j<ny; j++) {
      for (i=0; i<nx; i++) {
	x =  MapValue(iptr,i,j);
	MapValue(iptr,i,j) = pow(slope*(x-xmin), power);
      }
    }
    xmin = 0;
    xmax = 1;
  }
  
  fprintf(outstr,"P6\n%d %d\n255\n",nx,ny);
  for (j=ny-1; j>=0; j--) {
    for (i=0; i<nx; i++) {
      x = MapValue(iptr,i,j);
      get_rgb(x,xmin,xmax,&red,&green,&blue,Q8bit);
      fputc(red,outstr);
      fputc(green,outstr);
      fputc(blue,outstr);
    }
  }
  strclose(outstr);
}

int  get_lut(string fn, int n, real *r, real *g, real *b)
{
  real *coldat[3];
  int i, ncol, colnum[3];
  stream instr = stropen(fn,"r");

  colnum[0] = 1;    coldat[0] = r;
  colnum[1] = 2;    coldat[1] = g;
  colnum[2] = 3;    coldat[2] = b;
  ncol = get_atable(instr,3,colnum,coldat,MAXCOL);
  if (ncol < 0) 
    error("Problem %d reading %s; MAXCOL=%d",ncol,fn,MAXCOL);
  for (i=1; i<ncol; i++) {
    dRed[i-1]   = Red[i]   - Red[i-1];
    dGreen[i-1] = Green[i] - Green[i-1];
    dBlue[i-1]  = Blue[i]  - Blue[i-1];
  }
  return ncol;
}


void get_rgb(real data, double min, double max, int *R, int *G, int *B, bool Q8bit)
{
  int i;
  real x;
  
  if (min==max) {
    *R = *G = *B = (data > max ? 255 : 0);    
    return;
  }

  if (!Q8bit) {
    x = (data-min)*255.0/(max-min);
    if (x<=0.0 || x>=255.0) {         /* out of bounds */
      i=  (x <= 0.0 ? 0 : 255);
      *R = (int) (Red[i] * 255.0);
      *G = (int) (Green[i] * 255.0);
      *B = (int) (Blue[i] * 255.0);
      return;
    }
    i = (int) x;
    *R = (int)  ((Red[i]   + (x-i)*dRed[i])*255.0);
    *G = (int)  ((Green[i] + (x-i)*dGreen[i])*255.0);
    *B = (int)  ((Blue[i]  + (x-i)*dBlue[i])*255.0);
  } else {
    i = (int) ((data-min)*255.0/(max-min));
    x = i;
    if (i<0) i=0;
    if (i>255) i=255;
    *R = (int) (Red[i] * 255.0);
    *G = (int) (Green[i] * 255.0);
    *B = (int) (Blue[i] * 255.0);
  }
  dprintf(1,"%g [%g %g] -> (%g,%d) %d %d %d\n",data,min,max,x,i,*R, *G, *B);
  dprintf(2,"%g %g   %g %g    %g %g\n",Red[i],dRed[i],  Green[i],dGreen[i], Blue[i],dBlue[i]);
}

void stern_special(real data, real min, real max, int *pRed, int *pGreen, int *pBlue)
{
  real max_min = max - min;
  real a, b, c, d;
  real x, y;

  if(max_min > 0.0){
    x = (data - min)/(max_min);

    /* ========== Red ============ */
    a = 18.214286*x;
    b = 1.2857143 - 5.2040816*x;
    c = a < b ? a : b;
    y = c < 0.0 ? x : c;
    *pRed = (int)(256.0*y);
    *pRed = *pRed >   0 ? *pRed :   0;
    *pRed = *pRed < 255 ? *pRed : 255;

    /* ========== Green ========== */
    *pGreen = (int)(256.0*x);
    *pGreen = *pGreen >   0 ? *pGreen :   0;
    *pGreen = *pGreen < 255 ? *pGreen : 255;

    /* ========== Blue =========== */
    a = 1.9921875*x;
    b = 3.1333333 - 4.25*x;
    c = 3.8059701*x - 2.8059701;
    d = a < b ? a : b;
    y = d > c ? d : c;
    *pBlue = (int)(256.0*y);
    *pBlue = *pBlue >   0 ? *pBlue :   0;
    *pBlue = *pBlue < 255 ? *pBlue : 255;
  }
  else{
    *pRed = *pGreen = *pBlue = (data > max ? 255 : 0);
  }

  return;
}

