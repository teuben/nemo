/*
 * TABPEAK: find peaks - see also tablsqfit fit=peak
 *          
 *   28-may-2013   Created quick & dirty for ASTUTE               PJT
 *                        
 * 
 */

/**************** INCLUDE FILES ********************************/ 

#include <stdinc.h>
#include <strlib.h>
#include <getparam.h>
#include <moment.h>
#include <yapp.h>
#include <axis.h>
#include <mdarray.h>

/**************** COMMAND LINE PARAMETERS **********************/

string defv[] = {
    "in=???\n                     Input file name",
    "xcol=1\n			  X-Column",
    "ycol=2\n                     Y-column",
    "clip=0\n                     Only consider points above this",
    "nmax=100000\n                max size if a pipe",
    "VERSION=0.1\n		  28-may-2013 PJT",
    NULL
};

string usage = "peaks in a table";

string cvsid = "$Id$";

/**************** SOME GLOBAL VARIABLES ************************/

#ifndef MAXHIST
#define MAXHIST	1024
#endif

#ifndef MAXCOL
#define MAXCOL 256
#endif

#define MAXCOORD 16

local string input;			/* filename */
local stream instr;			/* input file */

local int col[2], ncol;

real *xcol, *ycol, *coldat[2];
local int    nmax;			/* lines to use at all */
local int    npt;			/* actual number of points */
real clip;


local void setparams(void);
local void read_data(void); 
local void peak_data(void);



/****************************** START OF PROGRAM **********************/

nemo_main()
{
    setparams();			/* read the parameters */
    read_data();
    peak_data();
}

local void setparams()
{
    input = getparam("in");
    col[0] = getiparam("xcol");
    col[1] = getiparam("ycol");
    clip = getrparam("clip");
    
    nmax = nemo_file_lines(input,getiparam("nmax"));
    if (nmax<1) error("Problem reading from %s",input);

    instr = stropen (input,"r");
}



local void read_data()
{
    int   i,j,k;
    
    dprintf(0,"Reading %d column(s)\n",ncol);
    xcol = (real *) allocate(sizeof(real)*nmax);
    ycol = (real *) allocate(sizeof(real)*nmax);
    coldat[0] = xcol;
    coldat[1] = ycol;
    ncol = 2;

    npt = get_atable(instr,ncol,col,coldat,nmax);        /* read it */
    if (npt == -nmax) {
    	warning("Could only read %d data",nmax);
    	npt = nmax;
    }
}


local void peak_data(void)
{
  int i,j;
  real mat[9], vec[3], sol[3], a[4];

  /* loop over all interior points and find peaks, fit local polynomial */

  for (i=1; i<npt-1; i++) {
    if (ycol[i] > clip && ycol[i]>ycol[i-1] && ycol[i]>ycol[i+1]) {
      lsq_zero(3,mat,vec);
      for (j=i-1; j<=i+1; j++) {
	a[0] = 1.0;
	a[1] = (xcol[j]-xcol[i]);
	a[2] = (xcol[j]-xcol[i]) * a[1];
	a[3] = ycol[j];
	lsq_accum(3,mat,vec,a,1.0);
      }
      lsq_solve(3,mat,vec,sol);
      dprintf(1,"Poly2 fit near i=%d (%g,%g)   %g %g %g\n",i+1,xcol[i],ycol[i],sol[0],sol[1],sol[2]);
      printf("%f %f \n",xcol[i] - sol[1]/(2*sol[2]),
	     sol[0]-sol[1]*sol[1]/(4*sol[2]));
    }
  }
}

