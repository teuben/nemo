/*
 * TABFFT: one-dim FFT's manipulator
 *
 *      3-apr-2008     1.0     written to test fftw/numrec  -      Peter Teuben
 */

/**************** INCLUDE FILES ********************************/ 
 
#include <stdinc.h>	
#include <getparam.h>

#include <fftw3.h>

#define MAXCOL  1

/**************** COMMAND LINE PARAMETERS **********************/

string defv[] = {                /* DEFAULT INPUT PARAMETERS */
    "in=???\n            Input file name (table)",
    "xcol=\n		 x coordinate column - REAL",
    "ycol=\n		 y coordinate column - IMAG",
    "nmax=0\n            max size of table to read",
    "duplicate=f\n       Double array to force complex conjugate",
    "VERSION=1.0\n	 3-apr-08 PJT",
    NULL
};

string usage = "tabular fft routines";

string cvsid="$Id$";


/**************** GLOBAL VARIABLES ************************/

local string input;				/* filename */

local int xcol[MAXCOL], ycol[MAXCOL], nxcol, nycol;	/* column numbers */

local real  *x[MAXCOL], *y[MAXCOL]; 			/* data from file */
local int    npt;				/* actual number of data points */

local int    nmax;				/* lines to allocate */

local bool Qdup;

void setparams();


/****************************** START OF PROGRAM **********************/

nemo_main()
{
    setparams();

    read_data();

    my_numrec();

}

void setparams()
{
    char *smin,*smax;
    int  i, j;
   
    input = getparam("in");             /* input table file */
    nxcol = nemoinpi(getparam("xcol"),xcol,MAXCOL);
    nycol = nemoinpi(getparam("ycol"),ycol,MAXCOL);
    if (nxcol < 1) error("Error parsing xcol=%s",getparam("xcol"));
    if (nycol < 0) error("Error parsing ycol=%s",getparam("ycol"));

    nmax = nemo_file_lines(input,getiparam("nmax"));
    dprintf(1,"Allocated %d lines for table\n",nmax);

    for (j=0; j<nxcol; j++)
        x[j] = (real *) allocate(sizeof(real) * (nmax+1));   /* X data */
    for (j=0; j<nycol; j++)
        y[j] = (real *) allocate(sizeof(real) * (nmax+1));   /* Y data */

    Qdup = getbparam("duplicate");

}

#define MVAL 		 64
#define MLINELEN	512

read_data()
{
    real *coldat[1+MAXCOL];
    int i, j, k, colnr[1+MAXCOL];
    stream instr;


    instr = stropen (input,"r");

		
    dprintf (2,"Reading datafile, xcol,ycol=%d..,%d,...\n",xcol[0],ycol[0]);
    for (j=0, k=0; j<nxcol; j++, k++) {
        colnr[k]  = xcol[j];
        coldat[k] = x[j];
    }
    for (j=0; j<nycol; j++, k++) {
        colnr[k]  = ycol[j];
        coldat[k] = y[j];
    }

    npt = get_atable(instr,nxcol+nycol,colnr,coldat,nmax);    /* get data */
    if (npt < 0) {
    	npt = -npt;
    	warning("Could only read first set of %d data",npt);
    }
    strclose(instr);
}


my_fftw()
{
  fftw_complex *in, *out;
  fftw_plan p;

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nmax);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nmax);
  p = fftw_plan_dft_1d(nmax, in, out, FFTW_FORWARD, FFTW_ESTIMATE);



  fftw_execute(p); /* repeat as needed */



  fftw_destroy_plan(p);
  fftw_free(in); 
  fftw_free(out);
}


#include "four1.c"

my_numrec()
{
  float *data;
  int i,j,n;

  n = nmax;

  data = (float *) allocate(nmax*4*sizeof(float));   /* 4 if Qdup,  2 if not */
  

  for (i=0, j=0; i<n; i++) {                     /* copy array for four1 */
    data[j++] = (nxcol ?  x[0][i] : 0.0);
    data[j++] = (nycol ?  y[0][i] : 0.0);
  }

  if (Qdup) {
    for (i=0, j=0; i<n; i++, j+=2) {             /* second complex conjugate part */
      data[4*n-2-j] =  data[j];                  /* real */
      data[4*n-1-j] = -data[j+1];                /* imag */
    }
    n *= 2;
  }

  for (i=0, j=0; i<n; i++, j+=2) {
    printf("SPEC1: %d %g %g\n",i, data[j], data[j+1]);
  }


  four1(data-1,n,1);

  for (i=0, j=0; i<n; i++, j+=2) {
    printf("LAG2: %d %g %g\n",i, data[j], data[j+1]);
  }

  four1(data-1,n,-1);

  for (i=0, j=0; i<n; i++, j+=2) {
    printf("SPEC3: %d %g %g\n",i, data[j]/n, data[j+1]/n);
  }


  
}

