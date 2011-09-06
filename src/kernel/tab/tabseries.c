/*
 * TABSERIES: discover repeated patterns
 *          
 *   2-sep-2011  Created quick & dirty for Cole Miller
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
    "xcol=1\n			  Column to use",
    "nw=3\n                       Coherence length to check for",
    "nmax=100000\n                max size if a pipe",
    "VERSION=0.1\n		  2-sep-2011 PJT",
    NULL
};

string usage = "discover repeated patterns";

string cvsid = "$Id$";

/**************** SOME GLOBAL VARIABLES ************************/



local string input;			/* filename */
local stream instr;			/* input file */

local int col;   			/* column number */
int   *w;
local int    nmax;			/* lines to use at all */
local int    npt;			/* actual number of points */
local int    nw;


local void setparams(void);
local void read_data(void); 
local void work_data(void);



/****************************** START OF PROGRAM **********************/

nemo_main()
{
  setparams();			/* read the parameters */
  read_data();
  work_data();
}

local void setparams()
{
  input = getparam("in");
  col = getiparam("xcol");
  nw = getiparam("nw");
  
  nmax = nemo_file_lines(input,getiparam("nmax"));
  if (nmax<1) error("Problem reading from %s",input);

  instr = stropen (input,"r");
}



local void read_data()
{
  int   i,j,k;
    
  w = (int *) allocate(sizeof(int)*nmax);
  npt = get_itable(instr,1,&col,&w,nmax);        /* read it */
  if (npt == -nmax) {
    warning("Could only read %d data",nmax);
    npt = nmax;
  }
}


local void work_data(void)
{
  int i,j,k;

  for (i=0; i<npt; i++) {
    for (j=i+nw; j<npt; j++) {
      for (k=0; k<npt; k++) {
	if (w[i+k] != w[j+k]) break;
      }
      if (k>nw) {
	printf("match: %d %d %d\n",i,j,k);
      }
    }
  }
}

