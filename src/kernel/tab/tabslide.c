/*
 * TABSLIDE:    slide through a large ascii table and do some work
 *
 *      11-jun-03   V0.1   template for Vanessa              PJT
 *
 * some performance numbers on a P4/1.6GHz
 * zcat 5body0.dat.gz | sum                                 10"
 * zcat 5body0.dat.gz | wc                                  56"
 * zcat 5body0.dat.gz | awk 'END{print NF,NR}'              35"
 * zcat 5body0.dat.gz | awk '{sum+=$1}END{print sum,NF,NR}' 40"
 * zcat 5body0.dat.gz | tabslide -                          34"
 * zcat 5body0.dat.gz | tabslide -  nmax=100000             35"
 */

#include <nemo.h>

string defv[] = {
    "in=???\n           input (table) file name (should be a pipe??)",
    "xcol=1\n           column(s) for x",
    "nmax=10000\n       Default max allocation",
    "VERSION=0.1\n      11-jun-03 PJT",
    NULL
};

string usage="template sliding through large ascii datasets";

#ifndef MAXCOL
#define MAXCOL 16
#endif

real *coldat[MAXCOL];
int colnr[MAXCOL];

stream instr, outstr;       /* input / output file */

int    nmax;                /* allocated space */
int    ntot;                /* total number of rows processed */
int    npt;                 /* number of rows in the buffer */
int    ncols;               /* number of colums */

void setparams(void);
int read_data(void);
void process_data(void);

real sum = 0.0;             /* some work related variables */


/****************************** START OF PROGRAM **********************/

nemo_main()
{
    setparams();
    while (read_data())
      process_data();
    dprintf(0,"\nProcessed %d lines\n",ntot);
    dprintf(0,"Sum = %g\n",sum);
}

/*
 *  get user interface and initialize variables 
 */

void setparams(void)
{
  int i;
  string inname = getparam("in");

  nmax = nemo_file_lines(inname,getiparam("nmax"));
  if (nmax<0) error("Error opening %s",inname);
  if (nmax==0) error("No data?");
  instr = stropen (inname,"r");
  
  ncols = nemoinpi(getparam("xcol"),colnr,MAXCOL);
  if (ncols<0) error("Illegal xcol= ncol=%d",ncols);
  
  for (i=0; i<ncols; i++)
    coldat[i] =  (real *) allocate(nmax * sizeof(real));
}

/* 
 * read a number of columns, and return number read, 0 for end of file
 */

int read_data(void)
{
  int i;

  npt = get_atable(instr,ncols,colnr,coldat,nmax);
  if (npt < 0) {
    npt = -npt;
  } 
  ntot += npt;
  dprintf(0,".",ntot);
  return npt;
}

/*
 * do some work, notice all variables are essentially global
 * not the worlds best style of programming 
 */

void process_data(void)
{
  int i, j;

  for (i=0; i<ncols; i++) {   /* loop over all col's */
    for (j=0; j<npt; j++) {   /* loop over all row's */
      sum += coldat[i][j];  /* just sum the data   */
      dprintf(1,"%d %d %g\n",i,j,coldat[i][j]);
    }
  }
}

