/*
 *  ASCII table routines
 *
 *	ASCII tables are columnar formatted, as long as first
 *	character in a line is not '#', '!' or ';', these are
 *	considered comment lines and automatically skipped
 *	as well as lines which consist of nothing but blank space
 *  
 *  get_Xtable:  return:   N > 0    N rows read, end of file
 *                         N = 0    nothing read, end of file
 *                         N < 0    abs(N) rows read, but no end of file yet 
 *
 *	 8-apr-92 PJT  Skip lines where #columns not correct
 *	22-jul-92 PJT  dprintf() comment lines 
 *	 6-aug-92 PJT  added get_ftable
 *	 3-feb-93 PJT  fixed bug in colfmt decl in get_ftable
 *	 4-mar-94 pjt  ansi headers
 *	15-jan-95 pjt  fixed SINGLEPREC problems
 *	24-jun-99 pjt  fixed reading tables important from DOS
 */

#include <stdinc.h>
#include <getparam.h>
#include <extstring.h>

#define MAXLINELEN  2048

extern string *burststring(string,string);

/*
 *  get_atable:  get table in memory, using free format
 *
 */

int get_atable(instr,ncol,colnr,coldat,ndat)
stream instr;                   /* in: input ascii file */
int ncol;                       /* in: number of columns to read */
int colnr[];                    /* in: column numbers to read */
real *coldat[];                 /* out: array of pointers to data */
int ndat;                       /* in: length of dat arrays */
{
    string *sp;
    int i, n, nr, nret, nline=0, npt=0;
    bool bad;
    char line[MAXLINELEN];
    real *dat;

    if (ndat<=0 || ncol<=0) error("Illegal ndat=%d ncol=%d",ndat,ncol);
    dprintf(2,"get_atable: parsing into %d columns\n",ncol);
    while (fgets(line,MAXLINELEN,instr) != NULL) {
        n = strlen(line);
        nline++;                        /* count number of lines read */
        if (line[n-1]=='\n') line[n-1]='\0';      /* patch line */
	if (line[0]=='#' || line[0]==';' || line[0]=='!') {
            dprintf(2,"%s\n",line);     
            continue;                       /* skip comment lines */
	}
        sp = burststring(line,", \t\r");     /* tokenize input line */
        n = xstrlen(sp,sizeof(string))-1;   /* number of items found */
        dprintf(3,"[%d] %s\n",n,line);
        if (n==0) continue;                 /* skip empty lines ? */
        if (npt >= ndat) { 
            npt = -ndat;                /* signal more data to come */
            break;
        }
        bad = FALSE;
        for (i=0; i<ncol; i++) {            /* process and fill each column */
            nr = colnr[i];		/* column number : >= 1 */
            dat = coldat[i];		/* pointer to data */
            if (nr>n) {		/* reference to an out of bounds column */
                warning("get_atable: skip line %d: %d columns, %d not present",
					nline,n,colnr[i]);
                bad = TRUE;
                break;
            }
            if (nr==0)			/* reference line number  */
                dat[npt] = nline;
            else if ( (nret=nemoinpr(sp[nr-1], &dat[npt], 1)) != 1) {
                warning("get_atable: line %d: error %d reading %s",
					nline,nret,sp[nr-1]);
                bad = TRUE;
                break;
            }
        } /* for (i) */
	freestrings(sp);
        if (bad) continue;
        npt++;                              /* count how much data filled */
    } /* while (fgets) */
    dprintf(1,"%d lines read from table file, %d data used \n",nline, npt);
    return npt;
}

/*
 *  get_ftable:  read fixed formatted table
 *
 */

int get_ftable(instr,ncol,colpos,colfmt,coldat,ndat)
stream instr;                   /* in: input ascii file */
int ncol;                       /* in: number of columns to read */
int colpos[];                   /* in: position to start reading */
string colfmt[];		/* in: format of that number */
real *coldat[];                 /* out: array of pointers to data */
int ndat;                       /* in: length of dat arrays */
{
    string *sp;
    int i, n, npos, nret, nline=0, npt=0, ecount=0;
    bool bad;
    char line[MAXLINELEN];
    real *dat;
    
    dprintf(2,"get_atable: parsing into %d columns\n",ncol);
    while (fgets(line,MAXLINELEN,instr) != NULL) {
        n = strlen(line);
        nline++;                        /* count number of lines read */
        if (line[n-1]=='\n') line[n-1]='\0';      /* patch line */
	if (line[0]=='#' || line[0]==';' || line[0]=='!') {
            dprintf(2,"%s\n",line);     
            continue;                       /* skip comment lines */
	}
        if (npt >= ndat) {                /* buffers full, signal */
            npt = -ndat;			/* that no eof */
            break;
        }
        bad = FALSE;
        for (i=0; i<ncol; i++) {            /* process and fill each column */
            npos = colpos[i];
            dat = coldat[i];
            if (npos == 0)
                dat[npt] = nline;
            else {
                if (sscanf(&line[npos-1],colfmt[i],&dat[npt]) != 1) 
                    bad = TRUE;
            }
        }
        if (bad) {
            ecount++;
            continue;
        }
        npt++;                              /* count how much data filled */
    } /* while (fgets) */
    if (ecount) warning("get_ftable: %d/%d/%d lines read were deemed bad",
        ecount, npt, nline);
    dprintf(1,"%d lines read from table file, %d data used \n",nline, npt);
    return npt;
}


#if defined(TOOLBOX)

#include <getparam.h>
#include <strlib.h>

string defv[] = {
    "in=???\n       Input table",
    "out=???\n      Output table",
    "colnr=\n       Input column numbers to extract, (free fmt) or:",
    "colpos=\n      Input column positions to extract (fixed fmt)",
    "colfmt=\n      Input format for fixed fmt case",
    "fmt=%g\n       Output format",
    "nmax=1000\n    Maximum number of rows (data) to read",
    "VERSION=1.0\n  6-aug-92",
    NULL,
};

string usage = "extract real table from a table";

#define MAXCOL   64

nemo_main()
{
    stream instr, outstr;
    int ndat;

    instr = stropen(getparam("in"),"r");
    ndat = getiparam("nmax");

    if (hasvalue("colnr")) {                /* free format */

        if (hasvalue("colpos")) warning("Ignoring colpos=");
        outstr = stropen(getparam("out"),"w");

        do_a(instr, outstr, ndat, getparam("colnr"),
             getparam("fmt"));

    } else if (hasvalue("colpos")) {        /* fixed format */

        if (hasvalue("colnr")) warning("Ignoring colnr=");
        outstr = stropen(getparam("out"),"w");

        do_f(instr, outstr, ndat, getparam("colpos"), getparam("colfmt"),
             getparam("fmt"));

    } else
        error("Neither colpos= nor colnr= has been supplied");
}


local string colfmt[MAXCOL];    /* output format */
local int    coli[MAXCOL];      /* fits either colnr or colpos */
local int    fmt[MAXCOL];       /* format for output */
local real  *coldat[MAXCOL];    /* pointers to data */

do_a(instr, outstr, ndat, scolnr, sfmt)
stream instr, outstr;
int ndat;
string scolnr, sfmt;
{   
    int ncol;

    ncol = nemoinpi(scolnr,coli,MAXCOL);
    
    my_alloc(ncol,ndat);
    error("free format not implemented");
}

do_f(instr, outstr, ndat, scolpos, scolfmt, sfmt)
stream instr, outstr;
int ndat;
string scolpos, scolfmt, sfmt;
{
    int ncol, n, i, j;
    string *sp, *burststring();
    real *dat;
    char fmt[30], *fmtftoc();

    ncol = nemoinpi(scolpos,coli,MAXCOL);       /* get column positions */

    sp = burststring(scolfmt," ,");             /* get input fmt's */
    if ( (n=xstrlen(sp,sizeof(string))) != ncol+1) 
        error("Not enough colfmt=; need %d, got %d",ncol,n);
    for (i=0; i<ncol; i++) {
        colfmt[i] = fmtftoc(sp[i]);             /* and convert/copy into colfmt */
    }

    sprintf(fmt,"%s ",sfmt);                    /* output format */

    my_alloc(ncol,ndat);                        /* allocate the data */
    do {                        /* loop while buffers full */
        n = get_ftable(instr,ncol,coli,colfmt,coldat,ndat);
        for (j=0; j<ABS(n); j++) {      /* output each row */
            for (i=0; i<ncol; i++) {        /* for each column */
                dat = coldat[i];
                fprintf(outstr,fmt,dat[j]);
            }
            fprintf(outstr,"\n");
        }
    } while (n<0);
    my_free(ncol,ndat);
}

my_alloc(ncol, ndat)
int ncol, ndat;
{
    int i;

    if (ncol > MAXCOL) error("my_alloc: Too many columns");

    for (i=0; i<ncol; i++)
        coldat[i] = (real *) allocate(ndat * sizeof(real));

}

my_free(ncol, ndat)
int ncol, ndat;
{
    int i;

    if (ncol > MAXCOL) error("my_free: Too many columns");

    for (i=0; i<ncol; i++)
        free ((char *)coldat[i]);
}

/* 
 *  convert fortran fmt statement to one in c, if possible
 *
 *  Fortran: (ANSI X3.9-1978 Section 13.2)
 *
 * Repeatable edit desciptors:      Non-repeatable:
 *  Iw                              'h1 h2 ... hn'
 *  Iw.m
 *  Fw.d
 *  Ew.d
 *  Ew.dEe
 *  Dw.d
 *  Gw.d
 *  Gw.dEe
 *  Lw
 *  A
 *  Aw
 *
 */

char *fmtftoc(s)
char *s;
{
    char fmt[32];

    if (s == NULL || *s == NULL) {  /* by default use %lf as format */
        strcpy(fmt,"%lf");
    } else if (*s == '%') {         /* C 'printf' format is passed 'as is' */
        sprintf(fmt,"%s",s);
    } else {                        /* assume fortran, try and parse */
        switch (*s) {
          case 'f':
          case 'F':
          case 'e':
          case 'E':
                s++;
                sprintf(fmt,"%%%slf",s);    /* also single precision */
                break;
          case 'g':
          case 'G':
          case 'd':
          case 'D':
                s++;
                sprintf(fmt,"%%%slf",s);    /* double precision */
                break;
          default:
                warning("Format specification \"%s\" not yet supported, assuming %lf",s);
                strcpy(fmt,"%lf");
                break;
        }
    }
    return scopy(fmt);    
}
#endif
