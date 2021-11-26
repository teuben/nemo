/*
 * MEANMED: NEMO version of a CARMA program with the same intent:
 *
 *      12-dec-2007   PJT    quick hack for running local 'quality'
 *
 */

#include <stdinc.h>	
#include <getparam.h>
#include <moment.h>



string defv[] = {              
    "infile=???\n        Input file name (table)",
    "maxpnt=10000\n      Max number of points that can be read",
    "carma=t\n           Special CARMA output",
    "VERSION=1.0\n	 12-dec-2007 PJT",
    NULL
};

string usage = "simple stats of all the numbers in a file";


local stream instr;


#ifndef MAX_LINELEN
#define MAX_LINELEN  2048
#endif

#ifndef MAXCOL
#define MAXCOL     1024
#endif

local char  line[MAX_LINELEN];
local real  data[MAXCOL];
local Moment m;



void nemo_main()
{
  int i,n;
  int maxpnt = getiparam("maxpnt");
  bool Qcarma = getbparam("carma");

  instr = stropen (getparam("infile"),"r");
  ini_moment(&m,2,maxpnt);

  for(;;) {
    if (fgets(line,MAX_LINELEN,instr) == NULL) 
      break;
    if (line[0] == '#') continue;
    n = strlen(line);
    if (line[n-1]=='\n') line[n-1]='\0';
    n = nemoinpr(line,data,MAXCOL);
    for (i=0; i<n; i++) {
      if (Qcarma && data[i]<=0.0) continue;
      accum_moment(&m,data[i],1.0);
    }
  }

  /* The CARMA version of meanmed outputs this:
   *   Npts       Median        Mean         Rms          Min          Max
   *    19        10.00         9.63         5.22         1.00        18.00
   */

  printf("  Npts       Median        Mean         Rms          Min          Max\n");
  printf("  %d ",n_moment(&m));
  printf("  %g ",median_moment(&m));
  printf("  %g ",mean_moment(&m));
  printf("  %g ",rms_moment(&m));
  printf("  %g ",min_moment(&m));
  printf("  %g ",max_moment(&m));
  printf("\n");


}

