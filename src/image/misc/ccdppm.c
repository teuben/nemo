/* 
 * CCDPPM:   create a PPM from a CCD image
 *
 *	16-dec-03  V1.0 created
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>
#include <image.h>
#include <table.h>

string defv[] = {
    "in=???\n       Input image filename",
    "out=???\n      Output PPM filename",
    "min=\n         Minimum overrride",
    "max=\n         Maximum overrride",
    "bad=\n         Use this as masking value to ignore data",
    "lut=\n         optional selection of a lut from NEMODAT/lut",
    "VERSION=1.0\n  16-dec-03 PJT",
    NULL,
};

string usage="convert image to a PPM file, with optional color cable";

string	infile, outfile;       		/* file names */
stream  instr, outstr;			/* file streams */

imageptr iptr=NULL;			/* will be allocated dynamically */
int    nx,ny,nz;			/* actual size of map */

real Red[256];
real Green[256];
real Blue[256];

void       get_rgb(real data, real min, real max, int *pR, int *pG, int *pB);
void      get_gray(real data, real min, real max, int *pR, int *pG, int *pB);
void stern_special(real data, real min, real max, int *pR, int *pG, int *pB);

int  get_lut(string, int, real*, real*, real*);

void nemo_main()
{
  int  i, j, red, green, blue, ncol;
  real x, xmin, xmax, bad, low;
  bool Qmin, Qmax, Qbad;
  
  instr = stropen (getparam("in"), "r");
  read_image (instr,&iptr);
  strclose(instr);
  outstr = stropen(getparam("out"), "w");
  
  nx = Nx(iptr);	
  ny = Ny(iptr);
  nz = Nz(iptr);
  Qmin = hasvalue("min");
  if (Qmin) xmin = getdparam("min");
  Qmax = hasvalue("max");
  if (Qmax) xmax = getdparam("max");
  Qbad = hasvalue("bad");
  if (Qbad) bad = getdparam("bad");

  if (hasvalue("lut"))
    ncol = get_lut(getparam("lut"),256,Red,Green,Blue);
  else {
    ncol = 256;
    for (i=0; i<ncol; i++)
      Red[i] = Green[i] = Blue[i] = i/255.0;
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
  }
  
  if (Qbad) {
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
  
  fprintf(outstr,"P6\n%d %d\n255\n",nx,ny);
  
  for (j=ny-1; j>=0; j--) {
    for (i=0; i<nx; i++) {
      x = MapValue(iptr,i,j);
#if 0
      get_rgb(x,xmin,xmax,&red,&green,&blue);
      get_grey(x,xmin,xmax,&red,&green,&blue);      
      stern_special(x,xmin,xmax,&red,&green,&blue);
#else
      get_rgb(x,xmin,xmax,&red,&green,&blue);
#endif
      fputc(red,outstr);
      fputc(green,outstr);
      fputc(blue,outstr);
    }
  }
}

void get_grey(real data, double min, double max, int *R, int *G, int *B)
{
  real x;
  if (min==max) {
    *R = *G = *B = (data > max ? 255 : 0);    
    return;
  }
  x = (data-min)*255.0/(max-min);
  if (x<=0.0 || x>=255.0) {         /* out of bounds */
    *R = *G = *B = (x <= 0.0 ? 0 : 255);
    return;
  }
  *R = *G = *B = (int) x;
  dprintf(0,"grey: %g\n",x);
}

void get_rgb(real data, double min, double max, int *R, int *G, int *B)
{
  int i;
  
  if (min==max) {
    *R = *G = *B = (data > max ? 255 : 0);    
    return;
  }
#if 0
  x = (data-min)*255.0/(max-min);
  if (x<=0.0 || x>=255.0) {         /* out of bounds */
    i=  (x <= 0.0 ? 0 : 255);
    *R = (int) (Red[i] * 255.0);
    *G = (int) (Green[i] * 255.0);
    *B = (int) (Blue[i] * 255.0);
    return;
  }
  i = (int) x;
  *R = (int)  (rgb[i*3+0] + (x-i)*der[i*3+0])*255.0;
  *G = (int)  (rgb[i*3+1] + (x-i)*der[i*3+1])*255.0;
  *B = (int)  (rgb[i*3+2] + (x-i)*der[i*3+2])*255.0;
#else
  i = (int) ((data-min)*255.0/(max-min));
  if (i<0) i=0;
  if (i>255) i=255;
  *R = (int) (Red[i] * 255.0);
  *G = (int) (Green[i] * 255.0);
  *B = (int) (Blue[i] * 255.0);
#endif
#if 0
  printf("%g [%g %g] -> %d %d %d %d\n",data,min,max,i,*R, *G, *B);
#endif
}

/* ========================================================================== */


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


/* ========================================================================== */

int  get_lut(string fn, int n, real *r, real *g, real *b)
{
  real *coldat[3];
  int ncol, colnum[3];
  stream instr = stropen(fn,"r");

  colnum[0] = 1;    coldat[0] = r;
  colnum[1] = 2;    coldat[1] = g;
  colnum[2] = 1;    coldat[2] = b;
  ncol = get_atable(instr,3,colnum,coldat,256);
  if (ncol < 0) error("Problem %d reading %s",ncol,fn);
  return ncol;
}
