/*
 * TABS: convert an ascii table to the SDAS/MIDAS-like tables
 *	 where all columns are separated by TABs
 *
 *      15-jul-94   created (see also: entab/detab)     pjt
 *
 */


#include <stdinc.h>
#include <getparam.h>
#include <ctype.h>

string defv[] = {
    "in=???\n           input file name(s)",
    "out=???\n          output file name",
    "col=???\n          column names, one for each, for header",
    "type=double\n      column types, one for each with last repeated",
    "delcol=\n          column numbers to skip at output (1=first)",
    "VERSION=1.0\n      16-jul-94 PJT",
    NULL,
};

string usage = "Table format converter";

extern string *burststring();



string input, output;			/* file names */
stream instr, outstr;			/* file streams */


#define MAXCOL          256             /* MAXIMUM number of columns */
#define MLINELEN       2048

bool   keepc[MAXCOL+1];                 /* columns to keep (t/f) */
string *colnames, *coltypes;
int    ncol, ntyp;

nemo_main()
{
    setparams();


    tabs();
}

setparams()
{
    int    delc[MAXCOL];                    /* columns to skip for output */
    int    i, ndelc;

    input = getparam("in");
    output = getparam("out");

    for (i=0; i<MAXCOL; i++)
        keepc[i]=TRUE;
    if (hasvalue("delcol")) {
        ndelc = nemoinpi(getparam("delcol"),delc,MAXCOL);
        if (ndelc<=0 || ndelc>MAXCOL)
            error("Too many columns given (%d) to delete (MAXCOL %d)",
			ndelc,MAXCOL);
        for (i=0; i<ndelc; i++)
            keepc[delc[i]-1] = FALSE;
    } 
    instr = stropen(input,"r");
    outstr = stropen (output,"w");

    colnames = burststring(getparam("col"),", \t");
    ncol = xstrlen(colnames,sizeof(string)) - 1;
    coltypes = burststring(getparam("type"),", \t");
    ntyp = xstrlen(coltypes,sizeof(string)) - 1;
    if (ncol<1) error("Need column names in: col=");
    if (ntyp<1) error("Need column types in: type=");

    for (i=0; i<ncol; i++) {
        if (i==0)
            fprintf(outstr,"%s",colnames[i]);
        else
            fprintf(outstr,"\t%s",colnames[i]);
    }
    fprintf(outstr,"\n");

    for (i=0; i<ntyp; i++) {
        if (i==0)
            fprintf(outstr,"%s",coltypes[i]);
        else
            fprintf(outstr,"\t%s",coltypes[i]);
    }
    for (i=ntyp; i<ncol; i++)
        fprintf(outstr,"\t%s",coltypes[ntyp-1]);
    fprintf(outstr,"\n");
    
}



tabs()
{
    char   line[MLINELEN];          /* input linelength */
    int    nval, i;
    string *outv;                   /* pointer to vector of strings to write */
    bool   tab;
        
    for(;;) {                       /* loop over all lines in file(s) */

        if (!get_line(instr, line)) return 0;
        dprintf(3,"LINE: (%s)\n",line);
        if(line[0] == '#') continue;		/* don't use comment lines */

        outv = burststring(line,"\t, ");
        nval = xstrlen(outv,sizeof(string)) - 1;
        if (nval != ncol) error("%d columns expected: \n%s\n",ncol,line);
	tab = FALSE;
        for (i=0; outv[i]; i++) {
            if (keepc[i]) {
                if (tab)  fputs("\t",outstr);
                fputs(outv[i],outstr);
                tab = TRUE;
            }
        }
        fputs("\n",outstr);
        freestrings(outv);
    } /* for(;;) */
}

