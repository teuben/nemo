/*
 *   FITSTAB:   output fits ASCII tables (XTENSION='TABLE   ')
 *
 *	8-apr-92    V1.0  Written       PJT
 *     15-apr-92    V1.1  added row=    PJT
 *     24-apr-92        a fixed bug in row= introduced earlier  PJT
 *     18-jan-92    V1.2  attempt to handle BINTABLE too  (fits.c)     PJT
 *     30-sep-96    V1.3  option to print out random group pars        pjt
 *
 *     10-nov-05    V1.5  add newline=
 * BUT PERHAPS WE SHOULD not boolean it, but with an integer give how often newlines are written
 *    
 */

#include <stdinc.h>
#include <getparam.h>
#include <fits.h>

string defv[] = {			/* Standard NEMO keyword+help */
    "in=???\n              Input fits table file",
    "hdu=0\n               Select which HDU file(s)? [0=all]",
    "select=data\n         Output mode {info,data,header,row,group}",
    "col=\n                Column names to select [all]",
    "row=\n                Row numbers to select [all]",
    "format=\n             Format for each row to use [TFORMnnn]",
    "maxrow=20000\n        Max row numbers to select",
    "random=f\n            Force random group?",
    "fnl=0\n               Frequency of newlines",
    "VERSION=1.5\n         10-nov-05 PJT",
    NULL,
};

string usage = "convert fits table or random groups to ascii table";

string cvsid="$Id$";

extern string *burststring(string,string);

void nemo_main()
{
    stream instr;
    int    i, n, nfile, maxrow, nrow, *row=NULL;
    string *col, select;
    struct fits_header fh;
    bool   scanopt();
    bool   Qrg = getbparam("random");
    int    fnl = getiparam("fnl");

    instr = stropen(getparam("in"),"r");    /* open input */
    
    if (hasvalue("col"))		/* select some columns */
        col = burststring(getparam("col"),", ");
    else				/* or select them all */
        col = NULL;
    select = getparam("select");
    nfile = getiparam("hdu");		/* need: nemoinpi() */
    maxrow = getiparam("maxrow");
    if (hasvalue("row")) {
        row = (int *) allocate((maxrow+1)*sizeof(int));
        nrow = nemoinpi(getparam("row"),row,maxrow);
        /* sorti(row,nrow); */
        if (nrow<0) {
            warning("Cannot use row selection - will use all");
            row = NULL;
        } else
            row[nrow] = 0;      /* signal end of list */
    }

    for (i=1;;i++) {			             /* try infinite loop */
       fts_zero(&fh);			             /* clean out header */
       n = fts_rhead(&fh,instr);	               /* read header */
       if (n<0)				              /* if no data (EOF) .. */
          break;			              /* ... quit */
       dprintf(0,"Working on FITS file %d\n",i);

       if ((nfile==0 || nfile==i)) {			/* if right file */
	 if (Qrg)
	   fts_pgroup(&fh,instr,col,select,row,fnl); 	/* try and print */
	 else
	   fts_ptable(&fh,instr,col,select,row,fnl); 	/* try and print */
       } else                                            /* Otherwise just */
	 fts_sdata(&fh,instr);		               /* skip the data */

       if (i==nfile)
            break;                                        /* done */
    }
    if (nfile>0 && i>=nfile && n<0)
        stop(1);                    /* see if early eof */

    strclose(instr);                    /* close files */
}
