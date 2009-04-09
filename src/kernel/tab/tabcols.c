/*
 * TABCOLS: select columns from a table
 *
 *      27-jan-00   created
 */

#include <stdinc.h>
#include <getparam.h>
#include <table.h>
#include <extstring.h>
#include <ctype.h>

string defv[] = {                /* DEFAULT INPUT PARAMETERS */
    "in=???\n           input file name(s)",
    "select=all\n       columns to select",
    "colsep=SP\n        Column separator (SP,TAB,NL)",    
    "out=-\n            output file name",
    "VERSION=2.0\n      9-apr-09 PJT",
    NULL
};

string usage = "Select column(s) from a table";


string input, output;			/* file names */
stream instr, outstr;			/* file streams */
int   ninput;				/* number of input files */

#ifndef MAX_COL
#define MAX_COL 256
#endif

#ifndef MAX_LINELEN
#define MAX_LINELEN  2048
#endif


#define MNEWDAT          80		/* space needed for one number */

int    keep[MAX_COL+1];                 /* column numbers to keep */
int    nkeep;                           /* actual number of skip columns */
int    maxcol = 0;                      /* largest column number */
char   colsep;                          /* character to separate columns */
bool   Qall;

local void setparams(void);
local void convert(stream , stream);
local void tab2space(char *);

extern  string *burststring(string, string);


nemo_main()
{
    setparams();
    instr  = stropen(input,"r");
    outstr = stropen (output,"w");
    convert (instr,outstr);
}

local void setparams(void)
{
    string select;                          /* which columns not to write */
    string separator;
    int i;

    input = getparam("in");
    output = getparam("out");

    select = getparam("select");
    if (select==NULL || *select=='\0' || streq(select,"all")) {
        Qall = TRUE;
        for (i=0; i<=MAX_COL; i++) keep[i] = i;
    } else {
        Qall = FALSE;
        nkeep = nemoinpi(select,keep,MAX_COL);
        if (nkeep<=0 || nkeep>MAX_COL)
            error("Too many columns given (%d) to keep (MAX_COL %d)",
      			nkeep,MAX_COL);
        for (i=0; i<nkeep; i++){
            if (keep[i]<0) error("Negative column number %d",keep[i]);
            maxcol = MAX(maxcol,keep[i]);
        }
    }

    separator = getparam("colsep");
    switch (*separator) {
        case 's':
        case 'S':       colsep =  ' ';   break;
        case 't':
        case 'T':       colsep = '\t';  break;
        case 'n':
        case 'N':       colsep = '\n';  break;
        case 'r':
        case 'R':       colsep = '\n';  break;
        case ':':       colsep =  ':';  break;
        case ',':       colsep =  ',';  break;
        default:        error("Illegal separator: s t n r : ,");
    }
}


local void convert(stream instr, stream outstr)
{
    char   line[MAX_LINELEN];          /* input linelength */
    int    i, nlines, noutv;
    string *outv;                   /* pointer to vector of strings to write */
    char   *cp, *seps=", \t";       /* column separators  */
        
    nlines=0;               /* count lines read so far */

    for (;;) {
        if (get_line(instr, line) < 0)      /* EOF */
            return; 					     

        dprintf(3,"LINE: (%s)\n",line);
        if (iscomment(line)) continue;
        
        nlines++;
        tab2space(line);                  /* work around a Gipsy (?) problem */

        outv = burststring(line,seps);
        noutv = xstrlen(outv,sizeof(string)) - 1;
        if (noutv < maxcol)
            error("Too few columns in input file (%d < %d)",noutv,maxcol);
        if (Qall) nkeep = noutv;
                        
        for (i=0; i<nkeep; i++) {
            if (keep[i] == 0)
                fprintf(outstr,"%d",nlines);
            else
                fprintf(outstr,"%s",outv[keep[i]-1]);
            if (i < nkeep-1) fprintf(outstr,"%c",colsep);
        }
        if (colsep != 'n') fprintf(outstr,"\n");    /* end of line */
    }
}
/*
 * small helper function, replaces tabs by spaces before processing.
 * this prevents me from diving into gipsy parsing routines and  fix
 * the problem there 
 * PJT - June 1998.
 */

local void tab2space(char *cp)
{
    while (*cp) {
        if (*cp == '\t') *cp = ' ';
        cp++;
    }
}
