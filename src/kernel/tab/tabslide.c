/*
 * TABSLIDE:    slide through a large dataset
 *
 *      11-jun-03   V0.1   template for Vanessa              PJT
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

#define MAXCOL 16

typedef struct column {
    int maxdat;     /* allocated length of data */          /* not used */
    int ndat;       /* actual length of data */             /* not used */
    real *dat;      /* pointer to data */
    int colnr;      /* column number this data came from */ /* not used */
} a_column;

int nxcol, nycol, xcolnr[MAXCOL];
a_column            xcol[MAXCOL];


stream instr, outstr;       /* input / output file */


int    nmax;                /* allocated space */
int    npt;                 /* actual number of points from table */


/****************************** START OF PROGRAM **********************/

nemo_main()
{
    setparams();
    read_data();
}


setparams()
{
    string inname = getparam("in");
    nmax = nemo_file_lines(inname,getiparam("nmax"));
    if (nmax<0) error("Error opening %s",inname);
    if (nmax==0) error("No data?");
    instr = stropen (inname,"r");

    nxcol = nemoinpi(getparam("xcol"),xcolnr,MAXCOL);
    if (nxcol<0) error("Illegal xcol= nxcol=%d",nxcol);
}


read_data()
{
    real *coldat[2*MAXCOL+1];
    int colnr[2*MAXCOL+1], ncols = 0, i, j, ntot=0;

    for (i=0; i<nxcol; i++) {
        coldat[ncols] = xcol[i].dat = (real *) allocate(nmax * sizeof(real));
        colnr[ncols] = xcolnr[i];        
        ncols++;
    }

    for (;;) {
      npt = get_atable(instr,ncols,colnr,coldat,nmax);
      if (npt < 0) {
        npt = -npt;
      } else if (npt==0) 
	break;
      ntot += npt;
      dprintf(0,"Processed %d\n",ntot);
    } 
	
}

