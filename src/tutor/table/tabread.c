/*
 * TABREAD: read a table, for the tutorial
 * 
 *      17-sep-2003     created                    PJT     
 *
 */

#include <nemo.h>


#define MAXCOL  256
#define MAXCOORD 16


string defv[] = {                /* DEFAULT INPUT PARAMETERS */
    "in=???\n            Input file name (table)",
    "xcol=1,2\n		 x coordinate column(s) to read",
    "nmax=100000\n       Hardcoded allocation space, if needed for piped data",
    "VERSION=1.0\n	 17-sep-03 PJT",
    NULL
};

string usage = "tutorial: example table reader";


nemo_main()
{
  real xmin, xmax;
  real *coldat[1+MAXCOL];
  int i, j, nxcol, colnr[1+MAXCOL], nmax, npt;
  string input;
  stream instr;

  input = getparam("in");                  /* input table filename */

  nmax = nemo_file_lines(input,getiparam("nmax"));
  dprintf(0,"Allocated %d lines for table\n",nmax);

  nxcol = nemoinpi(getparam("xcol"),colnr,MAXCOL);
  if (nxcol < 1) error("Error parsing xcol=%s",getparam("xcol"));

  dprintf(0,"Reading columns: ");
  for (j=0; j<nxcol; j++) {
    dprintf(0," %d",colnr[j]);
    coldat[j] = (real *) allocate(sizeof(real) * (nmax+1));
  }
  dprintf(0,"\n");

  npt = get_atable(instr,nxcol,colnr,coldat,nmax);    /* get data */
  if (npt < 0) {
    npt = -npt;
    warning("Could only read first set of %d data",npt);
  }
  dprintf(0,"Found %d lines in table\n",npt);
  
  /* find global min and max in all data */

  xmin = xmax =  coldat[0][0];
  for (i=0; i<npt; i++) {     
    for (j=0; j<nxcol; j++) {
      xmax=MAX(coldat[j][i],xmax);
      xmin=MIN(coldat[j][i],xmin);
    }
  }
}


