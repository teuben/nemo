/*   TABREAD: read an ASCII table in NEMO , for the tutorial
 * 
 *      17-sep-2003     created                    PJT     
 */

#include <nemo.h>

string defv[] = {
    "in=???\n         Input file name (table)",
    "xcol=1,2\n	      x coordinate column(s) to read (1=first)",
    "nmax=100000\n    Hardcoded allocation space, for piped data",
    "VERSION=1.0\n    19-sep-03 PJT",
    NULL
};

string usage = "tutorial: example table reader";

#define MAXCOL  16

nemo_main()
{
  int i, j, nxcol, colnr[MAXCOL], nmax, npt;
  real xmin, xmax, *coldat[MAXCOL];
  string input;
  stream instr;

  input = getparam("in");                             /* input table filename */
  nmax = nemo_file_lines(input,getiparam("nmax"));             /* count lines */
  dprintf(0,"Allocated %d lines for table\n",nmax);

  nxcol = nemoinpi(getparam("xcol"),colnr,MAXCOL);   /* colum numbers to read */
  if (nxcol < 1) error("Error parsing xcol=%s",getparam("xcol"));

  dprintf(0,"Reading %d column(s): ",nxcol);
  for (j=0; j<nxcol; j++) {          /* loop allocating space for each column */
    dprintf(0," %d",colnr[j]);
    coldat[j] = (real *) allocate(sizeof(real)*nmax);
  }
  dprintf(0,"\n");

  instr = stropen(input,"r");                                    /* open file */
  npt = get_atable(instr,nxcol,colnr,coldat,nmax);                /* get data */
  if (npt < 0) {
    npt = -npt;
    warning("Could only read first set of %d data",npt);
  }
  dprintf(0,"Found %d lines in table\n",npt);
  
  xmin = xmax =  coldat[0][0];          /* set the min/max to the first data */
  for (j=0; j<nxcol; j++) {                        /*  loop over all columns */
    for (i=0; i<npt; i++) {                      /* and loop over all points */
      xmax=MAX(coldat[j][i],xmax);
      xmin=MIN(coldat[j][i],xmin);
    }
  }
  dprintf(0,"MinMax in the data: %g %g\n",xmin,xmax);
}


