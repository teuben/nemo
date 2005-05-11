/* TABFILTER:   convolved a filter with a spectrum to get a flux
 *
 *
 *      10-may-05    first version, cloned off tabspline
 *
 */

#include <stdinc.h> 
#include <getparam.h>
#include <spline.h>

string defv[] = {
  "filter=???\n     Filter mnemonic (U,B,V,R,I,H,J,K) or full spectrum name",
  "spectrum=\n      Input spectrum table file (1=wavelength(A) 2=flux)",
  "tbb=\n           Black Body temperature, in case used",
  "xscale=1\n       Scale factor applied to input spectrum wavelength",
  "yscale=1\n       Scale factor applied to input spectrum flux",
  "step=1\n         Initial integration step (in Angstrom)",
  "VERSION=0.1\n    11-may-05 PJT",
  NULL,

};

string usage="flux derived from convolving a filter with a spectrum";

string cvsid="$Id$";


#define MAXZERO     64
#define MAXDATA  16384

/*
 * Planck curve,  output in ergs/cm2/s/A
 *
 * @index Planck curve, input: angstrom, kelvin   output: ergs/cm2/s/A
 */

real planck(real wavelen, real T_b)
{
  real w = wavelen/1e8;      /* w in cm now , for cgs units used here */
  real C1 = 3.74185e-5;
  real C2 = 1.43883;
  real x = C1/w/T_b;
  return C1/(w*sqr(sqr(w))*(exp(x)-1))*1e8;
}

nemo_main()
{
  int colnr[2];
  real *coldat[2], *xdat, *ydat, xmin, xmax, ymin, ymax;
  real fy[MAXZERO], fx[MAXZERO], xp[MAXDATA], yp[MAXDATA], x, y, xd, yd, dx;
  real tbb,sum;
  stream instr;
  string fmt, stype;
  char fmt1[100], fmt2[100];
  int i, j, n, nx, ny, nmax, izero;
  bool Qx, Qy, Qxall, Qyall;
  real *sdat;
  
  nmax = nemo_file_lines(getparam("filter"),MAXLINES);
  xdat = coldat[0] = (real *) allocate(nmax*sizeof(real));
  ydat = coldat[1] = (real *) allocate(nmax*sizeof(real));
  colnr[0] = 1;  /* wavelenght in angstrom */
  colnr[1] = 2;  /* normalized filter response [0,1] */
  
  instr = stropen(getparam("filter"),"r");
  n = get_atable(instr,2,colnr,coldat,nmax);
  
  for(i=0; i<n; i++) {
    dprintf(2,fmt,xdat[i],ydat[i]);
    if (i==0) {
      xmin = xmax = xdat[0];
      ymin = ymax = ydat[0];
    } else {
      xmax = MAX(xmax,xdat[i]);
      ymax = MAX(ymax,ydat[i]);
      xmin = MIN(xmin,xdat[i]);
      ymin = MIN(ymin,ydat[i]);
    }
  }
  dprintf(1,"Filter wavelength range: %g : %g\n",xmin,xmax);
  dprintf(1,"Filter response range: %g : %g\n",ymin,ymax);
  if (ydat[0]   != 0) warning("lower edge filter response not 0: %g",ydat[0]);
  if (ydat[n-1] != 0) warning("upper edge filter response not 0: %g",ydat[n-1]);
  dx = getdparam("step");
  if ((xmax-xmin)/100 < dx) {
    warning("Integration step %g in Angstrom too large, resetting to %g",
	    dx, (xmax-xmin)/100);
    dx = (xmax-xmin)/100;
  }
  
  /* setup a spline interpolation table into the filter */
  sdat = (real *) allocate(sizeof(real)*n*3);
  spline(sdat,xdat,ydat,n);
  
  if (hasvalue("tbb")) {                /* using a Planck curve */
    tbb = getdparam("tbb");
    
    sum = 0;
    for (x = xmin; x <= xmax; x += dx) {
      y = seval(x,xdat,ydat,sdat,n);    /* filter */
      dprintf(3,"%g %g\n",x,y);
      y *= planck(x,tbb);
      sum += y;
    }
    sum *= dx;
    printf("flux = %g   mag=%g\n",sum,-2.5*log10(sum));
  } else if (hasvalue("spectrum")) {
    warning("not implemented yet");
  } else
    warning("Either spectrum= or tbb= must be used");
}



