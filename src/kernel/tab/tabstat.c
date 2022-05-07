/*
 * TABSTAT: table statistics, like tabhist, but for multiple columns
 *          with iterative worst outlier removal -  no graphics output
 *
 *       9-dec-99   V1.0    Created
 *       6-jun-01   V1.1    renamed sigma= to nsigma= as in other programs PJT
 * 	 2-mar-01   V1.2    added sum to the output			   pjt
 *      24-jan-12   V1.4    also report min/max in sigma from mean         pjt
 *      16-jan-13   V1.5    added MAD                                      pjt
 *      10-oct-20   V1.7    ansi                                           pjt
 *      16-nov-21   V1.8    added qac= and robust=                         pjt
 *       1-dec-21   V1.9    with qac/robust keep the min/max from all data PJT
 *      23-apr-22   V2.0    new table V2 interface                         PJT
 *
 *  @todo:   xcol=0 should use the first data row to figure out all columns
 */

#include <stdinc.h>	
#include <getparam.h>
#include <moment.h>
#include <table.h>
#include <mdarray.h>

#define MAXCOL  256
#define MAXCOORD 16

string defv[] = {                /* DEFAULT INPUT PARAMETERS */
    "in=???\n            Input file name (table)",
    "xcol=1\n		 x coordinate column(s)",
    "iter=0\n            # iterations removing worst outlier",
    "nsigma=0\n          Multiple sigma above which to remove outliers",
    "verbose=t\n         Verbose output of all iterations?",
    "width=\n		 Width of columns to print",
    "median=t\n          Compute median too? (can be time consuming)",
    "mad=f\n             Compute MAD ?",
    "method=0\n          Method to remove outliers (0=fast 1=slow)",
#if 0    
    "nmax=100000\n       maximum number of data to be read if pipe",
#endif
    "xmin=\n             Set minimum ",
    "xmax=\n             Set maximum ",
    "bad=\n              Skip this bad value if one is given",
    "robust=f\n          robust stats?",
    "qac=f\n             QAC mode listing mean,rms,min,max",
    "label=\n            QAC label",
    "VERSION=2.2\n	 6-may-2022 PJT",
    NULL
};

string usage = "table column statistics";


local string input;				/* filename */
local stream instr;				/* input file */
local table  *tptr;

local int xcol[MAXCOL], nxcol; 		        /* column numbers */

local mdarray2  x;                              /*  x[col][row] */

local Moment m[MAXCOL];
local int    imaxdev[MAXCOL];
local int    npt;		                /* actual number of data points */
local int   *ix;

local bool   Qmedian;
local bool   Qverbose;
local bool   Qmad;
local bool   Qac;
local bool   Qrobust;
local bool   Qbad;
local real   badval;
local int    nmax;			 	 /* lines to allocate */
local int    width;
local char   outfmt[20];
local int    iter;
local real   fsigma;
local int    method;
local bool   Qmin, Qmax;
local real   xmin, xmax;
local string qac_label;

extern void sortptr (real *x ,int *idx, int n);

#define SWAP(a,b)   {real _t; _t=a; a=b; b=_t;}

void setparams(void);
void read_data(void);
void stat_data(void);
void out(string fmt);


void nemo_main(void)
{
    setparams();
    read_data();
    stat_data();
}

void setparams(void)
{
    int nrows, ncols;
   
    input = getparam("in");             /* input table file */
    instr = stropen (input,"r");
    tptr  = table_open(instr, 0);

    nrows = table_nrows(tptr);
    ncols = table_ncols(tptr);
    dprintf(0,"Table: %d x %d\n", nrows, ncols);

    nxcol = nemoinpi(getparam("xcol"),xcol,MAXCOL);
    if (nxcol == 0) {
      warning("selecting all columns");
      //nxcol = ncols;
      //if (ncols > MAXCOL) error("No room to select all (%d) columns; MAXCOL=%d", ncols,MAXCOL);
    } else if (nxcol < 1) {
      error("Error parsing xcol=%s   MAXCOL=%d",getparam("xcol"),MAXCOL);
    }
    nmax = nrows;

    Qverbose = getbparam("verbose");
    Qmedian = getbparam("median");
    Qmad = getbparam("mad");
    Qac = getbparam("qac");
    Qrobust = getbparam("robust");
    if (hasvalue("width")) {
        width = getiparam("width");
        sprintf(outfmt,"%%%ds",width);
    } else
        sprintf(outfmt,"%%s");
    Qmin = hasvalue("xmin");
    Qmax = hasvalue("xmax");
    Qbad = hasvalue("bad");
    if (Qmin) xmin = getrparam("xmin");
    if (Qmax) xmax = getrparam("xmax");
    if (Qbad) badval = getrparam("bad");
    dprintf(1,"format=%s\n",outfmt);
    fsigma = getrparam("nsigma");
    iter = getiparam("iter");
    method = getiparam("method");
    if (hasvalue("label"))
      qac_label = getparam("label");
    else
      qac_label = getparam("in");
}

void read_data(void)
{
    // x[col][row] will be the data
    x = table_md2cr(tptr, nxcol, xcol, 0, 0);
    npt = tptr->nr;
    if (nxcol == 0) nxcol = tptr->nc;
}


void stat_data(void)
{
    int i, j, ndat, imax, kmin, kmax;
    real median, mean, sigma, d, dmax, rrange[2];
    char fmt[20];
    
    ix = (int *) allocate(sizeof(int)*npt);     /* pointer array */

    ndat = 0;
    if (Qmad || Qac || Qrobust) ndat = npt;

    for (j=0; j<nxcol; j++) {           /* initialize moments for all data */
        ini_moment(&m[j],4,ndat);
        for (i=0; i<npt; i++) {                          /* loop over rows */
	  if (Qbad && x[j][i]==badval) continue;
	  if (Qmin && x[j][i]<xmin) continue;
	  if (Qmax && x[j][i]>xmax) continue;
	  accum_moment(&m[j],x[j][i],1.0);
        }
    }

    // simpler one line output
    if (Qac) {   
      for (j=0; j<nxcol; j++) {
	if (Qrobust) {
	  compute_robust_moment(&m[j]);
	  robust_range(&m[j], rrange);
	  printf("QAC_STATS: %s %g %g %g %g  %g %g  %d\n",
		 qac_label, mean_robust_moment(&m[j]), sigma_robust_moment(&m[j]),
		 min_moment(&m[j]), max_moment(&m[j]),		 
		 sum_moment(&m[j]), sratio_moment(&m[j]), n_robust_moment(&m[j]));
	} else
	  printf("QAC_STATS: %s %g %g %g %g  %g %g  %d\n",
		 qac_label, mean_moment(&m[j]), sigma_moment(&m[j]),
		 min_moment(&m[j]), max_moment(&m[j]),
		 sum_moment(&m[j]), sratio_moment(&m[j]),n_moment(&m[j]));
      }
      return;
    }

    do {                                /* iteration loop to reject outliers */

        if (Qverbose || iter==0) {          /* always print out last iter */
                                            /* and in verbose all iters */
            printf("npt:    ");
            for (j=0; j<nxcol; j++) {
                sprintf(fmt," %d",n_moment(&m[j]));
                out(fmt);
            }
            printf("\n");
    
            printf("min:    ");
            for (j=0; j<nxcol; j++) {
                sprintf(fmt," %g",min_moment(&m[j]));
                out(fmt);
            }
            printf("\n");

            printf("max:    ");
            for (j=0; j<nxcol; j++) {
                sprintf(fmt," %g",max_moment(&m[j]));   
                out(fmt);
            }
            printf("\n");
    
            printf("sum:    ");
            for (j=0; j<nxcol; j++) {
  	        sprintf(fmt," %g",sum_moment(&m[j]));
                out(fmt);
            }   
            printf("\n");

            printf("mean:   ");
            for (j=0; j<nxcol; j++) {
                sprintf(fmt," %g",mean_moment(&m[j]));
                out(fmt);
            }   
            printf("\n");

            printf("disp:   ");
            for (j=0; j<nxcol; j++) {
                sprintf(fmt," %g",sigma_moment(&m[j]));        
                out(fmt);
            }
            printf("\n");
	    if (Qmad) {
	      printf("mad:    ");
	      for (j=0; j<nxcol; j++) {
                sprintf(fmt," %g",mad_moment(&m[j]));        
                out(fmt);
	      }
	      printf("\n");
	    }

            printf("skew:   ");
            for (j=0; j<nxcol; j++) {
                sprintf(fmt," %g",skewness_moment(&m[j]));        
                out(fmt);
            }
            printf("\n");

            printf("kurt:  ");
            for (j=0; j<nxcol; j++) {
                sprintf(fmt," %g",kurtosis_moment(&m[j]));        
                out(fmt);
            }
            printf("\n");


	    printf("min/sig:");
            for (j=0; j<nxcol; j++) {
	      d = (min_moment(&m[j]) - mean_moment(&m[j]))/sigma_moment(&m[j]);
	      sprintf(fmt," %g",d);
	      out(fmt);
            }
            printf("\n");

	    printf("max/sig:");
            for (j=0; j<nxcol; j++) {
	      d = (max_moment(&m[j]) - mean_moment(&m[j]))/sigma_moment(&m[j]);
	      sprintf(fmt," %g",d);
	      out(fmt);
            }
            printf("\n");


            if (Qmedian) {
                printf("median: ");
                for (j=0; j<nxcol; j++) {
                    sortptr(x[j],ix,npt);
                    kmin = 0;
                    kmax = npt-1;
                    if ((kmax-kmin+1)%2) {
                        median = x[j][ix[kmin+(kmax-kmin)/2]];
                    } else {
                        median = 0.5 * (x[j][ix[kmin+(kmax-kmin+1)/2]] +
                                        x[j][ix[kmin+(kmax-kmin+1)/2-1]]);
                    }
                    sprintf(fmt," %g",median);
                    out(fmt);
                }
                printf("\n");
            }
        }


        if (iter) {    /* if another iteration needed to reject worst outlier */
            if (Qverbose || iter==0) printf("\niter %d:",iter);
            for (j=0; j<nxcol; j++) {
                mean = mean_moment(&m[j]);
                sigma = sigma_moment(&m[j]);
                dmax = 0.0;
                imax = npt-1;
                if (method==0) {
                  for (i=0; i<npt; i++) {      /* fast method */
                    d = (x[j][i]-mean)/sigma;
                    d = ABS(d);
                    if (d > dmax) {
                        imax = i;
                        dmax = d;
                    }
                  }
                } else if (method==1) {
                  for (i=0; i<npt; i++) {      /* slow method : only works if all points in sample */
                    decr_moment(&m[j],x[j][i],1.0);
                    d = (x[j][i]-mean_moment(&m[j]))/
                            sigma_moment(&m[j]);
                    d = ABS(d);
                    if (d > dmax) {
                        imax = i;
                        dmax = d;
                    }
                    accum_moment(&m[j],x[j][i],1.0);
                  }
                } else
                  error("Illegal outlier removal method %d",method);
                dprintf(1,"Swapping %g and %g in %d and %d\n",
                        x[j][imax],x[j][npt-1],imax,npt-1);
                SWAP(x[j][imax],x[j][npt-1]);
                imaxdev[j] = imax;                  /* remember where */
                sprintf(fmt," %g",dmax);
                if (Qverbose || iter==0) out(fmt);
            }
            if (Qverbose || iter==0) {
                printf("\n");
                printf("idx  %d:",iter);
                for (j=0; j<nxcol; j++) {
                    sprintf(fmt," %d",imaxdev[j]+1);
                    out(fmt);
                }
                printf("\n");
            }
            npt--;

            dprintf(1,"Redoing %d columns %d rows\n",nxcol,npt);
            for (j=0; j<nxcol; j++) {       /* redo, for min/max */
                reset_moment(&m[j]);
                for (i=0; i<npt; i++) {
                    accum_moment(&m[j],x[j][i],1.0);
                }
            }
        }
    } while (iter--);
    free(ix);
}


void out(string fmt)
{
    printf(outfmt,fmt);   
}
