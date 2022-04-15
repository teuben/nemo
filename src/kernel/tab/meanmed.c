/*
 * MEANMED: NEMO version of a CARMA program with the same intent:
 *
 *      12-dec-2007   PJT    1.0 quick hack for running local 'quality'
 *       4-mar-2022   PJT    1.1 converted to table V2 (initial steps) + cleanup
 *
 */

#include <stdinc.h>	
#include <getparam.h>
#include <moment.h>
#include <table.h>


string defv[] = {              
    "infile=???\n        Input file name (table)",
    "maxpnt=10000\n      Max number of points that can be read",
    "carma=t\n           Special CARMA output",
    "VERSION=1.1\n	 4-mar-2022 PJT",
    NULL
};

string usage = "simple stats of all the numbers in a file";


#ifndef MAXCOL
#define MAXCOL     1024
#endif


void nemo_main()
{
  bool Qcarma = getbparam("carma");
  int maxpnt = getiparam("maxpnt");
  real  data[MAXCOL];
  stream instr;
  tableptr tp;
  Moment m;
  string s;
  int i,k=0,n, mode = 0;

  instr = stropen(getparam("infile"),"r");
  tp = table_open(instr, mode);

  ini_moment(&m,2,maxpnt);

  while(1) {
    if (!(s=table_line0(tp)))        // @todo   table_line0 
      break;
    if (s[0] == '#') continue;
    n = strlen(s);
    if (s[n-1]=='\n') s[n-1]='\0';
    n = nemoinpr(s,data,MAXCOL);
    for (i=0; i<n; i++) {
      if (Qcarma && data[i]<=0.0) continue;
      k++;
      if (k==maxpnt) error("too many points, increase maxpnt=%d",maxpnt);
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

